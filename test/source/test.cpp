#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>
// #include <spread/spatidatamanager.h>
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
    // CSpreadAnalyse sp;
    // sp.elevationPath = path;
    // sp.InitEnvironment();
    CFieldStrengthAnalyse *analyse = new CFreeSpaceAnalyse;
    // 测试坐标
    // 103.556378  31.566910
    Station station;
    station.lontitude = 103.556378;
    station.latitude = 31.566910;
    station.frequency = 107.2;
    // CStations stations;
    // stations.AddStation(&station);
    analyse->elevationPath = path;
    analyse->needComputeAll = true;
    analyse->pStations->AddStation(&station);
    bool flag = false;
    static_cast<CFreeSpaceAnalyse *>(analyse)->FieldStrengthAnalyse(
        "/home/shigp/data", RasterCreateFileType::rcftTiff, &flag);
    std::cout << "return flag is:" << flag << std::endl;
    IFileFloatArray *arr = analyse->pData;
    long cout = 0;
    arr->GetSize(&cout);
    for (long i = 0; i < cout; i++) {
      float *p = arr->GetRasterDataArray();
      std::cout << p[i] << std::endl;
    }
  }
}
