#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>
#include <pcldemo/pcldemo.h>
#include <pcldemo/version.h>

#include <string>
TEST_CASE("PCLDemo") {
  using namespace pcldemo;

  PCLDemo pcldemo("Tests");
  std::string aa = pcldemo.program_version();
  std::cout << aa << std::endl;
  MESSAGE(aa);
}