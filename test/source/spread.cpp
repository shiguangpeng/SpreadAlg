#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>
#include <spread/spread.h>
#include <spread/version.h>

#include <string>
TEST_CASE("Spread") {
  using namespace spread;

  Spread spread("Tests");
  std::string aa = spread.program_version();
  std::cout << aa << std::endl;
  MESSAGE(aa);
}