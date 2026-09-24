#pragma once

// Dependency-free unit test framework: TEST registers a case, CHECK* record failures and keep
// going, REQUIRE* abort the current case. Run the binary with a substring to filter cases.

#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace samaya::test {

struct Case {
  const char* name;
  std::function<void()> body;
};

inline std::vector<Case>& registry() {
  static std::vector<Case> cases;
  return cases;
}

inline int& failures() {
  static int count = 0;
  return count;
}

struct Registrar {
  Registrar(const char* name, std::function<void()> body) {
    registry().push_back({name, std::move(body)});
  }
};

struct AbortCase {};

inline void report(const char* file, int line, const std::string& message) {
  ++failures();
  std::cerr << file << ":" << line << ": FAILED " << message << "\n";
}

}  // namespace samaya::test

#define SAMAYA_TEST_CONCAT2(a, b) a##b
#define SAMAYA_TEST_CONCAT(a, b) SAMAYA_TEST_CONCAT2(a, b)

#define TEST(name)                                                              \
  static void SAMAYA_TEST_CONCAT(test_fn_, name)();                             \
  static const ::samaya::test::Registrar SAMAYA_TEST_CONCAT(test_reg_, name)(   \
      #name, &SAMAYA_TEST_CONCAT(test_fn_, name));                              \
  static void SAMAYA_TEST_CONCAT(test_fn_, name)()

#define CHECK(cond)                                                             \
  do {                                                                          \
    if (!(cond)) ::samaya::test::report(__FILE__, __LINE__, "CHECK(" #cond ")"); \
  } while (0)

#define CHECK_EQ(a, b)                                                          \
  do {                                                                          \
    const auto& va_ = (a);                                                      \
    const auto& vb_ = (b);                                                      \
    if (!(va_ == vb_)) {                                                        \
      std::ostringstream os_;                                                   \
      os_ << "CHECK_EQ(" #a ", " #b "): " << va_ << " != " << vb_;              \
      ::samaya::test::report(__FILE__, __LINE__, os_.str());                    \
    }                                                                           \
  } while (0)

#define CHECK_NEAR(a, b, tol)                                                   \
  do {                                                                          \
    const double va_ = (a);                                                     \
    const double vb_ = (b);                                                     \
    if (!(std::fabs(va_ - vb_) <= (tol))) {                                     \
      std::ostringstream os_;                                                   \
      os_ << "CHECK_NEAR(" #a ", " #b "): " << va_ << " vs " << vb_;            \
      ::samaya::test::report(__FILE__, __LINE__, os_.str());                    \
    }                                                                           \
  } while (0)

#define CHECK_THROWS(expr, exception_type)                                      \
  do {                                                                          \
    bool thrown_ = false;                                                       \
    try {                                                                       \
      (void)(expr);                                                             \
    } catch (const exception_type&) {                                           \
      thrown_ = true;                                                           \
    }                                                                           \
    if (!thrown_) {                                                             \
      ::samaya::test::report(__FILE__, __LINE__,                                \
                             "CHECK_THROWS(" #expr ", " #exception_type ")");   \
    }                                                                           \
  } while (0)

#define REQUIRE(cond)                                                           \
  do {                                                                          \
    if (!(cond)) {                                                              \
      ::samaya::test::report(__FILE__, __LINE__, "REQUIRE(" #cond ")");         \
      throw ::samaya::test::AbortCase{};                                        \
    }                                                                           \
  } while (0)
