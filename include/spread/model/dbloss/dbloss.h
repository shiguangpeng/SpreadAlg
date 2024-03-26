/**
 * @copyright ALL COPYRIGH RESERVED BY SHIGP
 * @author shigp
 */

#ifndef INCLUDE_SPREAD_MODEL_DBLOSS_DBLOSS_H_
#define INCLUDE_SPREAD_MODEL_DBLOSS_DBLOSS_H_

#include <gdal_priv.h>
#include <spatdata/datatype.h>
#include <spread/model/combine.h>
#include <spread/model/spreadanalyse.h>
#include <spread/spread.h>

using datatype::Station;
using spatdata::IAnalyseEnvironment;

// 前置声明，避免编译报错
struct ICombineAnalyse;

namespace dbloss {
// DBLossElement结构体，父结构体，保存绕射传播的各种参数
struct IDBLossElement {
 public:
  // 是否复杂模型
  virtual void getComplicatedModel(bool *isComplicated) = 0;
  // 设置CombineAnalyse
  virtual void setCombineAnalyse(ICombineAnalyse *combineAnalyse) = 0;
  virtual void prepareAnalyseEnvi(void) = 0;
  virtual double GetDBLoss(Station station, OGRPoint point) = 0;
  virtual double GetDBLossRev(Station station, OGRPoint point) = 0;
};

// 实现类声明
class CDBLossElement {
 public:
  CDBLossElement()
      : Cols(0),
        Rows(0),
        CellSize(0),
        XMin(0.0),
        YMax(0.0),
        NoData(0),
        pEnvi(nullptr),
        pElevData(nullptr),
        pAnalyse(nullptr),
        hm(0) {}

  virtual ~CDBLossElement() {
    // 清理
    if (pEnvi) {
      delete pEnvi;
      pEnvi = nullptr;
    }
    if (pElevData) {
      delete pElevData;
      pElevData = nullptr;
    }
    if (pAnalyse) {
      delete pAnalyse;
      pAnalyse = nullptr;
    }
  }

 protected:
  int64_t Cols;
  int64_t Rows;
  double CellSize;
  double XMin;
  double YMax;
  double NoData;
  IAnalyseEnvironment *pEnvi;
  IFileFloatArray *pElevData;
  ICombineAnalyse *pAnalyse;
  float hm;

 protected:
  /**
   * @brief 准备分析环境
   */
  void PrepareAnalyseEnvi();
  /**
   * @brief 设置联合分析环境
   * @param newVal
   */
  void SetCombineAnalyse(ICombineAnalyse *newVal);
  /**
   * @brief 获得高程值
   * @param X 横坐标
   * @param Y 纵坐标
   * @return 获取传入的坐标处的栅格值
   */
  float GetZValue(double X, double Y);
};

}  // namespace dbloss

#endif  // INCLUDE_SPREAD_MODEL_DBLOSS_DBLOSS_H_
