// samaya command-line interface.
//
//   samaya [options] model.mps
//
// With --json, a single JSON object summarizing the run is printed as the last line of stdout;
// bench/harness.py relies on that.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <string_view>

#include "samaya.hpp"

namespace {

constexpr int kExitOk = 0;
constexpr int kExitReadError = 1;
constexpr int kExitUsage = 2;

void print_usage(std::FILE* out) {
  std::fprintf(out,
               "Usage: samaya [options] <model.mps>\n"
               "\n"
               "Options:\n"
               "  --stats               Print model statistics and exit without solving\n"
               "  --json                Print a JSON run summary as the last line of output\n"
               "  --time-limit <sec>    Wall-clock time limit\n"
               "  --threads <n>         Worker threads (0 = all cores)\n"
               "  --mip-gap <gap>       Relative MIP gap tolerance\n"
               "  --lp-method <m>       auto | dual | primal | barrier | pdlp | concurrent\n"
               "  --no-presolve         Disable presolve\n"
               "  --gpu                 Use GPU kernels when available\n"
               "  --log-level <0-3>     Verbosity (default 1)\n"
               "  --version             Print the version and exit\n"
               "  --help                Print this message and exit\n");
}

bool parse_lp_method(std::string_view s, samaya::LpMethod& out) {
  using samaya::LpMethod;
  if (s == "auto") out = LpMethod::kAuto;
  else if (s == "dual") out = LpMethod::kDualSimplex;
  else if (s == "primal") out = LpMethod::kPrimalSimplex;
  else if (s == "barrier") out = LpMethod::kBarrier;
  else if (s == "pdlp") out = LpMethod::kPdlp;
  else if (s == "concurrent") out = LpMethod::kConcurrent;
  else return false;
  return true;
}

std::string json_string(std::string_view s) {
  std::string out = "\"";
  for (const char c : s) {
    switch (c) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\t': out += "\\t"; break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) {
          char buf[8];
          std::snprintf(buf, sizeof buf, "\\u%04x", c);
          out += buf;
        } else {
          out += c;
        }
    }
  }
  return out + "\"";
}

std::string json_number(double v) {
  if (!std::isfinite(v)) return "null";
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.17g", v);
  return buf;
}

}  // namespace

int main(int argc, char** argv) {
  samaya::Params params;
  bool stats_only = false;
  bool json = false;
  std::string path;

  for (int i = 1; i < argc; ++i) {
    const std::string_view arg = argv[i];
    const auto next_value = [&](const char* name) -> const char* {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "samaya: %s requires a value\n", name);
        std::exit(kExitUsage);
      }
      return argv[++i];
    };
    const auto next_number = [&](const char* name) -> double {
      const char* text = next_value(name);
      char* end = nullptr;
      const double v = std::strtod(text, &end);
      if (end == text || *end != '\0') {
        std::fprintf(stderr, "samaya: invalid value '%s' for %s\n", text, name);
        std::exit(kExitUsage);
      }
      return v;
    };

    if (arg == "--help" || arg == "-h") {
      print_usage(stdout);
      return kExitOk;
    } else if (arg == "--version") {
      std::printf("samaya %s\n", samaya::version());
      return kExitOk;
    } else if (arg == "--stats") {
      stats_only = true;
    } else if (arg == "--json") {
      json = true;
    } else if (arg == "--time-limit") {
      params.time_limit = next_number("--time-limit");
    } else if (arg == "--threads") {
      params.threads = static_cast<int>(next_number("--threads"));
    } else if (arg == "--mip-gap") {
      params.mip_rel_gap = next_number("--mip-gap");
    } else if (arg == "--lp-method") {
      const char* method = next_value("--lp-method");
      if (!parse_lp_method(method, params.lp_method)) {
        std::fprintf(stderr, "samaya: unknown LP method '%s'\n", method);
        return kExitUsage;
      }
    } else if (arg == "--no-presolve") {
      params.presolve = false;
    } else if (arg == "--gpu") {
      params.use_gpu = true;
    } else if (arg == "--log-level") {
      params.log_level = static_cast<int>(next_number("--log-level"));
    } else if (!arg.empty() && arg.front() == '-') {
      std::fprintf(stderr, "samaya: unknown option '%s'\n", argv[i]);
      print_usage(stderr);
      return kExitUsage;
    } else if (path.empty()) {
      path = arg;
    } else {
      std::fprintf(stderr, "samaya: more than one model file given\n");
      return kExitUsage;
    }
  }
  if (path.empty()) {
    print_usage(stderr);
    return kExitUsage;
  }

  const auto read_start = std::chrono::steady_clock::now();
  samaya::Model model;
  try {
    model = samaya::read_mps(path);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "samaya: error reading '%s': %s\n", path.c_str(), e.what());
    if (json) {
      std::printf("{\"instance\":%s,\"status\":\"read_error\",\"message\":%s}\n",
                  json_string(path).c_str(), json_string(e.what()).c_str());
    }
    return kExitReadError;
  }
  const double read_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - read_start).count();

  const samaya::ModelStats stats = samaya::compute_stats(model);
  if (params.log_level > 0 || stats_only) {
    samaya::print_stats(std::cout, model, stats);
    std::cout << "Read time " << read_seconds << " s\n";
    std::cout.flush();
  }

  samaya::Result result;
  if (!stats_only) result = samaya::Solver(params).solve(model);

  if (json) {
    std::printf(
        "{\"instance\":%s,\"name\":%s,\"class\":\"%s\",\"rows\":%d,\"cols\":%d,\"nnz\":%lld,"
        "\"integers\":%d,\"status\":\"%s\",\"objective\":%s,\"dual_bound\":%s,"
        "\"read_seconds\":%s,\"solve_seconds\":%s,\"nodes\":%lld,\"message\":%s}\n",
        json_string(path).c_str(), json_string(model.name).c_str(),
        samaya::to_string(model.problem_class()), stats.rows, stats.cols,
        static_cast<long long>(stats.nnz), stats.integers, samaya::to_string(result.status),
        json_number(result.objective).c_str(), json_number(result.dual_bound).c_str(),
        json_number(read_seconds).c_str(), json_number(result.solve_seconds).c_str(),
        result.nodes, json_string(result.message).c_str());
  }
  return kExitOk;
}
