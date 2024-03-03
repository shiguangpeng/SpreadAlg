#include "spread/spatidatamanager.h"

#include <gdal_priv.h>

#include <iostream>

using namespace spatidatamanager;

SpatialDataManagement::SpatialDataManagement() {}
SpatialDataManagement::~SpatialDataManagement() {}

RasterDataSource::RasterDataSource() {
  // 注册所有驱动
  GDALAllRegister();
}
RasterDataSource::~RasterDataSource() {}

GDALDataset* RasterDataSource::OpenRaster(const char* path) {
  // 打开模式，只读

  return static_cast<GDALDataset*>(GDALOpen(path, GA_ReadOnly));
}
