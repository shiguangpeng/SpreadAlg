/**
 * @copyright ALL COPYRIGH RESERVED BY SHIGP
 * @author shigp
 */

#ifndef INCLUDE_SPREAD_MODEL_DBLOSS_DBLOSS_H_
#define INCLUDE_SPREAD_MODEL_DBLOSS_DBLOSS_H_

#include <gdal_priv.h>
#include <spatdata/datatype.h>
#include <spread.h>
// #include <spread/model/combine.h>
#include <spread/spreadbase.h>

using datatype::Station;
using spatdata::IAnalyseEnvironment;
using spread::ICombine;

namespace dbloss {
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
  ICombine *pAnalyse;
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
  void SetCombineAnalyse(ICombine *newVal);
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
