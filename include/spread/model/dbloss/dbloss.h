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

// 声明有这个类
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
  virtual double_t GetDBLoss(Station station, OGRPoint point) = 0;
  virtual double_t GetDBLossRev(Station station, OGRPoint point) = 0;
};
}  // namespace dbloss

#endif  // INCLUDE_SPREAD_MODEL_DBLOSS_DBLOSS_H_
