#include <gdal_priv.h>
#include "spread/spatidatamanager.h"


using namespace spatidatamanager;

#pragma region SpatialDataManagement
SpatialDataManagement::SpatialDataManagement(){}

SpatialDataManagement::~SpatialDataManagement(){}
#pragma endregion

#pragma region RasterDataSource
RasterDataSource::RasterDataSource(){}
RasterDataSource::~RasterDataSource(){}

int RasterDataSource::OpenRaster(const char *path, GDALDataset* poDataset){
    GDALAllRegister();
    const char * filePath = path;
    // 打开模式，只读
    GDALAccess access = GA_ReadOnly;
    poDataset = GDALDataset::FromHandle(GDALOpen(filePath, access));
    if (poDataset== nullptr)
    {
        return -1;
    }
    return 0;
}
#pragma endregion


