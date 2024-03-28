// All Copyright reserved. shigp
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>
#include <spatdata/spatdata.h>
#include <spread.h>
#include <spread/model/fieldstrength.h>
#include <spread/model/freespace/freespace.h>
#include <spread/spreadbase.h>

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
    using datatype::RasterCreateFileType;
    using datatype::Station;
    using spatdata::CGDALRasterReaderByPixel;
    using spatdata::CStations;
    using spatdata::IFileFloatArray;
    using spread::spreadbase::combine::freespace::CFreeSpace;
    using spread::spreadbase::fieldstrength::CFieldStrength;

    std::string path = "/home/shigp/data/test_prj.tif";
    CFieldStrength *analyse = new CFreeSpace;
    // 测试坐标
    Station station;
    station.x = 11502076.828;
    station.y = 3505543.756;
    station.frequency = 107.2;
    // station.stationHeight = 10;

    station.stationHeight = 10;
    // station.stationHeight = 100;
    station.power = 10 * log10(200);
    station.stationHeight = 30;
    station.power = 10 * log10(100);  // 一般不使用W，而使用dBw, 故有此转化
    station.freqThreshold = 40;

    analyse->elevationPath = path;
    analyse->needComputeAll = true;
    analyse->pStations = new CStations;
    analyse->pStations->AddStation(&station);
    bool flag = false;
    static_cast<CFreeSpace *>(analyse)->FieldStrengthAnalyse(
        "/home/shigp/data", RasterCreateFileType::rcftTiff, &flag);
    std::cout << "return flag is:" << flag << std::endl;
    IFileFloatArray *arr = analyse->pData;
    int64_t count = 0;
    arr->GetSize(&count);

    float *data = arr->GetRasterDataArray();
    // // 打开文件进行写操作
    // std::ofstream outfile;
    // std::string savePath = "/home/shigp/data/result_90mheight_test.csv";
    // outfile.open(savePath, std::ios::out);
    // // 写入数组内容到文件
    // for (int64_t i = 0; i < count; i++) {
    //   if (i == count - 1) {
    //     outfile << data[i] << std::endl;
    //   } else {
    //     outfile << data[i] << ",";
    //   }
    // }

    // 测试，gdal保存栅格数据
    GDALAllRegister();
    GDALDriver *driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    GDALDataset *saveData =
        driver->Create("/home/shigp/data/out_offset.tif", analyse->cols,
                       analyse->rows, 1, GDALDataType::GDT_Float32, nullptr);

    GDALDataset *openFileData = static_cast<GDALDataset *>(
        GDALOpen(analyse->elevationPath.c_str(), GA_ReadOnly));
    double geoTransform[6];
    openFileData->GetGeoTransform(geoTransform);
    saveData->SetGeoTransform(geoTransform);

    const OGRSpatialReference *ref = openFileData->GetSpatialRef();
    saveData->SetSpatialRef(ref);

    for (u_int32_t i = 0; i < (u_int32_t)sizeof(data) / sizeof(*data); i++) {
      if (data[i] < 0) {
        data[i] = 0;
      }
    }
    saveData->GetRasterBand(1)->SetNoDataValue(0);

    // GDALRasterBand *band = saveData->GetRasterBand(1);
    // CPLErr err =
    //     band->RasterIO(GF_Write, 0, 0, analyse->cols, analyse->rows, data,
    //                    analyse->cols, analyse->rows, GDT_Float32, 0, 0);

    CPLErr err = saveData->RasterIO(
        GF_Write, 0, 0, analyse->cols, analyse->rows, data, analyse->cols,
        analyse->rows, GDT_Float32, 1, nullptr, 0, 0, 0, nullptr);
    if (err != CE_None) {
      std::cout << "栅格写入失败!"
                << "错误代码：" << err << std::endl;
    }

    saveData->Close();
    GDALDestroy();

    // 关闭文件
    // outfile.close();
    // 清理内存
    delete analyse;
  }
}
