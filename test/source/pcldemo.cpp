/**
 * @copyright all copyright reserved
 * @author shigp
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>
#include <pcldemo.h>
#include <pcldemo/version.h>

#include <string>
TEST_CASE("PCLDemo") {
  using pcldemo::PCLDemo;

  PCLDemo pcldemo("Tests");
}