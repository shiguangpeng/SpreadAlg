#include "spread/spatidatamanager.h"

#include <gdal_priv.h>

#include <iostream>

using namespace spatidatamanager;

#pragma region SpatialDataManagement

SpatialDataManagement::SpatialDataManagement() {}

SpatialDataManagement::~SpatialDataManagement() {}
#pragma endregion

#pragma region "RasterDataSource定义"
RasterDataSource::RasterDataSource() {}
RasterDataSource::~RasterDataSource() {}

GDALDataset* RasterDataSource::OpenRaster(const char* path) {
  // 打开模式，只读
  GDALAllRegister();
  return static_cast<GDALDataset*>(GDALOpen(path, GA_ReadOnly));
}
#pragma endregion