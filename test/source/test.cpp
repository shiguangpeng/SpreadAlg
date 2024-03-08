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

  TEST_CASE("自由传播模型调用测试") {
    using namespace spread;
    std::string path = "/home/shigp/data/test.tif";
    CSpreadAnalyse sp = CSpreadAnalyse();
    sp.elevationPath = path;
    sp.InitEnvironment();

    CFreeSpaceAnalyse analyse;
    bool flag;
    analyse.FieldStrengthAnalyse("/home/shigp/data", RasterCreateFileType::rcftTiff, &flag);
    std::cout << "result is:" << flag << std::endl;
  }
}
