#include <fmt/format.h>
#include <spread/spread.h>

#include <iostream>

#include "spread/spatidatamanager.h"

using namespace spread;
using namespace spatidatamanager;

CSpreadAnalyse::CSpreadAnalyse() {}
CSpreadAnalyse::~CSpreadAnalyse() {}

bool CSpreadAnalyse::InitEnvironment() {
  // 通过打开的栅格数据，初始化类相关成员。
  const char* path = elevationPath.c_str();
  RasterDataSource ptr = RasterDataSource();
  GDALDataset* data = ptr.OpenRaster(path);
  if (data == nullptr) {
    errorInfo = "栅格数据打开失败。";
    return false;
  }
  // 使用GDAL数据指针初始化栅格属性
  // 行列数
  cols = data->GetRasterXSize();
  rows = data->GetRasterYSize();
  int bandNum = data->GetRasterCount();
  // 波段从1开始
  // for (int i = 1; i <= bandNum; i++) {
  //   GDALRasterBand* band = data->GetRasterBand(i);
  //   int x = 0;
  //   int y = 0;
  // }

  // 结束前关闭数据集对象
  GDALClose(data);
  GDALDestroy();
  errorInfo = "栅格数据打开成功。";
  return true;
}
