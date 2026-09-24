#pragma once

#include <chrono>
#include <cstdarg>

namespace samaya {

// Minimal leveled logger writing to stdout. Level semantics match Params::log_level.
class Logger {
 public:
  explicit Logger(int level) : level_(level) {}

  int level() const { return level_; }
  bool enabled(int level) const { return level <= level_; }

  // printf-style; a trailing newline is added.
  void log(int level, const char* fmt, ...) const
#if defined(__GNUC__)
      __attribute__((format(printf, 3, 4)))
#endif
      ;

 private:
  int level_;
};

class Timer {
 public:
  Timer() : start_(std::chrono::steady_clock::now()) {}
  double seconds() const {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
  }

 private:
  std::chrono::steady_clock::time_point start_;
};

}  // namespace samaya
