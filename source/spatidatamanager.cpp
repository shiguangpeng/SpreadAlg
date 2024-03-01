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

int RasterDataSource::OpenRaster(const char* path, GDALDataset* ptr) {
  // 打开模式，只读
  GDALAllRegister();
  ptr = static_cast<GDALDataset*>(GDALOpen(path, GA_ReadOnly));
  char** metadata = ptr->GetMetadata();
  // 打开tiff文件

  std::cout << "metadata: " << metadata << std::endl;
  if (ptr == NULL) {
    return -1;
  }
  return 0;
}
#pragma endregion