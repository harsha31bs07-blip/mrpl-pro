#pragma once

#include <cstdint>

namespace samaya {

enum class Status : std::uint8_t {
  kNotSolved,
  kOptimal,
  kInfeasible,
  kUnbounded,
  kInfeasibleOrUnbounded,
  kTimeLimit,
  kIterationLimit,
  kNodeLimit,
  kNumericalError,
  kInvalidModel,
  kNotImplemented,
};

const char* to_string(Status status);

}  // namespace samaya
