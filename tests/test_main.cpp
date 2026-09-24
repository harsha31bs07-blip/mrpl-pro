#include <cstring>
#include <exception>
#include <iostream>

#include "test_framework.hpp"

int main(int argc, char** argv) {
  using namespace samaya::test;
  const char* filter = argc > 1 ? argv[1] : nullptr;
  int ran = 0;
  int failed_cases = 0;
  for (const Case& c : registry()) {
    if (filter != nullptr && std::strstr(c.name, filter) == nullptr) continue;
    const int before = failures();
    try {
      c.body();
    } catch (const AbortCase&) {
    } catch (const std::exception& e) {
      report(c.name, 0, std::string("unexpected exception: ") + e.what());
    }
    ++ran;
    if (failures() != before) {
      ++failed_cases;
      std::cerr << "[FAIL] " << c.name << "\n";
    } else {
      std::cout << "[ OK ] " << c.name << "\n";
    }
  }
  std::cout << ran - failed_cases << "/" << ran << " test cases passed\n";
  return failed_cases == 0 && ran > 0 ? 0 : 1;
}
