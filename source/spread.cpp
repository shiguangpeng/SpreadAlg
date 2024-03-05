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
  // 仿射变换参数
  double_t geoInfo[6];
  data->GetGeoTransform(geoInfo);
  xMin = geoInfo[0];
  yMax = geoInfo[3];
  cellSize = geoInfo[1];
  // 获取当前tif图的所有波段数量
   int bandNum = data->GetRasterCount();
  // 波段从1开始，dem这种结果图层也只有一个图层
    int *dataTypeFlag = nullptr;
    double_t bandNoData = 0.0;
  for (int i = 1; i <= bandNum; i++) {
    GDALRasterBand* band = data->GetRasterBand(i);
    GDALDataType dataType = band->GetRasterDataType();
    // 根据GetNoDataValue()方法的说明，特定数据类型需要使用其他方法
    if (dataType == GDALDataType::GDT_Int64)
    {
      bandNoData = band-> GetNoDataValueAsInt64(dataTypeFlag);
    }else if(dataType == GDALDataType::GDT_UInt64){
      bandNoData = band->GetNoDataValueAsUInt64(dataTypeFlag);
    }else {
      bandNoData = band->GetNoDataValue(dataTypeFlag);
    }
    // todo: 就项目来说，像dem这种结果tiff图，土地利用类型等，一般只有一个波段。这里先break;
    break;
  }
  if (dataTypeFlag)
  {
    noData = bandNoData;
  }
  // 结束前关闭数据集对象
  GDALClose(data);
  GDALDestroy();
  errorInfo = "栅格数据打开成功。";
  return true;
}
