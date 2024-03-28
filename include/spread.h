/**
 * @copyright ALL COPYRIGH RESERVED BY SHIGP
 * @author shigp
 * @note all relevant of spread interface defined here.
 */

#ifndef INCLUDE_SPREAD_H_
#define INCLUDE_SPREAD_H_

#include <error.h>
#include <spatdata/datatype.h>
#include <spatdata/spatdata.h>
#include <spread/spreadbase.h>

#include <iostream>
#include <string>
#include <vector>

#include "/usr/local/include/gdal_priv.h"

using std::string;
using std::vector;

using datatype::RasterCreateFileType;
using spatdata::IStations;
using spread::spreadbase::AnalyseEnvironment;

namespace spread {
// 在顶层命名空间中使用前置声明
struct ICombine;
// ISpread、IFieldStrength、ICombine、IFreeSpace
struct ISpreadBase {
 public:
  /**
   * @brief 获取分析环境
   */
  virtual void getAnalyseEnvi(IAnalyseEnvironment **pVal) = 0;

  /**
   * @brief 设置分析环境
   */
  virtual void put_AnalyseEnvi(IAnalyseEnvironment *newVal) = 0;

  /**
   * @brief 获取dem路径
   */
  virtual void get_ElevationPath(string *path) = 0;

  /**
   * @brief 设置dem路径
   */
  virtual void put_ElevationPath(string path) = 0;

  /**
   * @brief 获取错误信息
   */
  virtual void get_ErrorInfo() = 0;
};

/**
 * @brief IFieldStrength接口，继承了ISpreadBase
 */
struct IFieldStrength : public ISpreadBase {
  /**
   * @brief 获取要分析的台站列表
   * @param stations
   */
  virtual void get_Stations(IStations **stations) = 0;

  /**
   * @brief 设置要分析的台站列表
   * @param stations
   */
  virtual void put_Stations(IStations *stations) = 0;

  /**
   * @brief 获取hm参数
   */
  virtual void get_hm(float *hm) = 0;

  /**
   * @brief 设置hm参数
   */
  virtual void put_hm(float hm) = 0;

  /**
   * @brief 获取改正的传播强度
   */
  virtual void get_OffsetDB(float *db) = 0;

  /**
   * @brief 设置改正的传播强度
   */
  virtual void put_OffsetDB(float db) = 0;

  /**
   * @brief 场强分析
   * @param outPutPath
   * @param type
   * @param pVal
   */
  virtual void FieldStrengthAnalyse(string outPutPath,
                                    RasterCreateFileType type, bool *pVal) = 0;

  virtual void FieldStrengthCounting(string outPutPath,
                                     RasterCreateFileType type, bool *pVal) = 0;

  /**
   * @brief 获取要分析的子区域
   * @param pVal
   */
  virtual void get_SubExtent(string *pVal) = 0;

  /**
   * @brief 设置要分析的子区域
   * @param newVal
   */
  virtual void put_SubExtent(string newVal) = 0;

  /**
   * @brief 单点场强分析
   * @param lon 经度
   * @param lat 纬度
   * @param pVal
   */
  virtual void FieldStrengthAnalyseByPoint(double lon, double lat,
                                           vector<double> **pVal) = 0;

  /**
   * @brief 清理资源
   */
  virtual void DestroyData() = 0;
};

/**
 * DBLossElement结构体，父结构体，保存绕射传播的各种参数
 */
struct IDBLossElement {
 public:
  // 是否复杂模型
  virtual void getComplicatedModel(bool *isComplicated) = 0;
  // 设置CombineAnalyse
  virtual void setCombineAnalyse(ICombine *combineAnalyse) = 0;
  virtual void prepareAnalyseEnvi(void) = 0;
  virtual double GetDBLoss(Station station, OGRPoint point) = 0;
  virtual double GetDBLossRev(Station station, OGRPoint point) = 0;
};

/**
 * @brief ICombine接口，继承了IFieldStrength
 * @note spread相关的继承体系，没有规范到对应的实现类，从ICombine开始
 * 由CCombine继承，项目中所有的ICCombine基类类型均要指向其子类的实现类。
 */
struct ICombine : IFieldStrength {
  virtual ~ICombine();

  /**
   * @brief 添加其他影像数据
   * @param path
   * @param progressName
   */
  virtual void AddOtherRaster(string path, string progressName) = 0;

  /**
   * @brief 获取其他影像数据
   * @param pVal
   */
  virtual void GetOtherRaster(IFileFloatArray **pVal) = 0;

  /**
   * @brief 获取高程数据
   * @param pVal
   */
  virtual void get_ElevationData(IFileFloatArray **pVal) = 0;

  /**
   * @brief 获取文件读取高程信息对象
   * @param pVal
   */
  virtual void get_Elevation(IGDALRasterReaderByFileArray **pVal) = 0;

  /**
   * @brief 获取遥感影像中的nodata的值
   * @param pVal
   */
  virtual void get_NoData(double *pVal) = 0;

  /**
   * @brief 添加绕射造成的强度损失
   * @param el
   */
  virtual void AddDBLossElement(IDBLossElement *el) = 0;
};
}  // namespace spread
#endif  // INCLUDE_SPREAD_H_
