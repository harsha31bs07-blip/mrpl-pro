#include "core/log.hpp"

#include <cstdio>

namespace samaya {

void Logger::log(int level, const char* fmt, ...) const {
  if (!enabled(level)) return;
  va_list args;
  va_start(args, fmt);
  std::vfprintf(stdout, fmt, args);
  va_end(args);
  std::fputc('\n', stdout);
  std::fflush(stdout);
}

}  // namespace samaya
