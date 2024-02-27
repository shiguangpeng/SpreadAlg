#pragma once

#include <string>
#include <spread/version.h>

namespace spread {
  /**
   * @brief A basic class for spread model 
   * 
   */
  class Spread {
    std::string name;

  public:
    /**
     * @brief Creates a new spread
     * @param name the name to greet
     */
    Spread(std::string name);

    /**
     * @brief indicate the version of program
     * @return a string containing the version info
     */
    std::string program_version() const;
    const char* test_gdal();
  };

  /**
   * @brief 场强分析算法的基类
  */
  class DEMAnalyse{
    public:
      DEMAnalyse();
    virtual ~DEMAnalyse();
    protected:
    bool init_environment();
    protected:
    struct analyse_env_st *p_envi;
    std::string elevation_path;





  };

  /**
   * @brief 数据源管理基类，用于对数据的操作，如：栅格、矢量以及空间数据库，该类以及子类都是对GDAL能力的封装
   * @attention GDAL >= 3.5(ver) supported
  */
  class DataSourceManagement {
    public:
     DataSourceManagement();
  };

  class RasterManagement: public DataSourceManagement {

  };

  class VectorManagement: public DataSourceManagement{

  };

}  // namespace spread
