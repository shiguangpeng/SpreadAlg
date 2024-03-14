// All Copyright reserved. shigp
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>
#include <spread/spread.h>

#include <fstream>
#include <iostream>
#include <string>

/// @brief 传播算法测试SUITE
TEST_SUITE("传播算法测试用例") {
  // TEST_CASE("SPAT") {
  //   using namespace spatidatamanager;
  //   GDALDataset* result =
  //   RasterDataSource::OpenRaster("F:\\work\\spreadmodel\\test.tif");
  // }

  TEST_CASE("自由传播模型调用测试") {
    using spread::CFieldStrengthAnalyse;
    using spread::CFreeSpaceAnalyse;
    using spread::CStations;
    using spread::IFileFloatArray;
    using spread::RasterCreateFileType;
    using spread::Station;

    std::string path = "/home/shigp/data/test_prj.tif";
    CFieldStrengthAnalyse *analyse = new CFreeSpaceAnalyse;
    // 测试坐标
    Station station;
    station.x = 11502076.828;
    station.y = 3505543.756;
    station.frequency = 107.2;
    station.stationHeight = 30;
    station.power = 10 * log10(100);
    station.freqThreshold = 40;
    station.dem = 499;
    analyse->elevationPath = path;
    analyse->needComputeAll = true;
    analyse->pStations = new CStations();
    analyse->pStations->AddStation(&station);
    bool flag = false;
    static_cast<CFreeSpaceAnalyse *>(analyse)->FieldStrengthAnalyse(
        "/home/shigp/data", RasterCreateFileType::rcftTiff, &flag);
    std::cout << "return flag is:" << flag << std::endl;
    IFileFloatArray *arr = analyse->pData;
    int64_t count = 0;
    arr->GetSize(&count);

    float *data = arr->GetRasterDataArray();
    // 打开文件进行写操作
    std::ofstream outfile;
    outfile.open("/home/shigp/data/result.csv", std::ios::out);

    // 写入数组内容到文件
    for (int64_t i = 0; i < count; i++) {
      if (i == count - 1) {
        outfile << data[i] << std::endl;
      } else {
        outfile << data[i] << ",";
      }
    }
    // 关闭文件
    outfile.close();
    // 清理内存
    delete analyse;
  }
}
