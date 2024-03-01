#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>
// #include <spread/spread.h>
#include <spread/spatidatamanager.h>
#include <spread/version.h>

#include <string>
#include <gdal_priv.h>

// TEST_CASE("Spread") {
//   using namespace spread;

//   Spread spread("Tests");
//   std::string aa = spread.program_version();
//   const char *p = spread.test_gdal();
//   MESSAGE(*p);
// }

/// @brief 空间数据管理的测试用例
TEST_SUITE("空间数据管理的测试用例"){
  TEST_CASE("SPAT"){
    using namespace spatidatamanager;

    RasterDataSource source = RasterDataSource();
    GDALDataset* point = NULL;
    int flag = source.OpenRaster("F:\\work\\spreadmodel\\test.tif", point);
    MESSAGE("打开结果：{}", flag);
  }
}
