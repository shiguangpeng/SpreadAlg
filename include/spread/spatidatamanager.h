#ifndef SPREADALG_INCLUDE_SPATIAL_SUPPORT_
#define SPREADALG_INCLUDE_SPATIAL_SUPPORT_

#include <error.h>

#include "gdal_priv.h"

namespace spatidatamanager {
  /**
   * @brief
   * 数据源管理基类，用于对数据的操作，如：栅格、矢量以及空间数据库，该类以及子类都是对GDAL能力的封装
   * @attention GDAL >= 3.5(ver) supported project
   */
  class SpatialDataManagement {
  public:
    /**
     * @brief 初始化时初始化所有支持的矢量和栅格数据格式的驱动
     * @see https://gdal.org/drivers/raster/index.html and
     * https://gdal.org/drivers/vector/index.html
     */
    SpatialDataManagement();
    ~SpatialDataManagement();
  };

  class RasterDataSource : public SpatialDataManagement {
  public:
    RasterDataSource();
    ~RasterDataSource();

    /// @brief 根据传入的路径，打开对应的栅格数据，返回数据的GDALDataset指针
    /// @param filePath
    /// @return GDALDataset* 该文件的GDALDataset指针
    GDALDataset* OpenRaster(const char* filePath);
  };

  /**
   * @brief 封装了GDAL矢量数据的相关操作
   */
  class VectorDataSoure : public SpatialDataManagement {};

}  // namespace spatidatamanager
#endif