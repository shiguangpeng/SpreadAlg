#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>
#include <spread/spatidatamanager.h>
#include <spread/spread.h>

#include <string>

/// @brief 空间数据管理的测试用例
TEST_SUITE("空间数据管理的测试用例") {
  // TEST_CASE("SPAT") {
  //   using namespace spatidatamanager;
  //   GDALDataset* result = RasterDataSource::OpenRaster("F:\\work\\spreadmodel\\test.tif");
  // }

  TEST_CASE("InitEnvironment测试") {
    using namespace spread;
    std::string path = "F:\\work\\spreadmodel\\test.tif"; // "/home/shigp/Downloads/test.tif"; // F:\work\spreadmodel
    CSpreadAnalyse sp = CSpreadAnalyse();
    sp.elevationPath = path;
    sp.InitEnvironment();
  }
}
