/**
 * @copyright ALL COPYRIGH RESERVED BY SHIGP
 * @author shigp
 */

#ifndef INCLUDE_SPREAD_MODEL_FIELDSTRENGTHANALYSE_H_
#define INCLUDE_SPREAD_MODEL_FIELDSTRENGTHANALYSE_H_

#include <spatdata/datatype.h>
#include <spatdata/spatdata.h>
#include <spread/model/spreadanalyse.h>

#include <string>
#include <vector>

using datatype::RasterCreateFileType;
using datatype::Station;
using spatdata::IAnalyseEnvironment;
using spatdata::IStations;
using spread::spreadanalyse::CSpreadAnalyse;
using std::string;
using std::vector;

namespace spread {
namespace spreadanalyse {
namespace fieldstrengthanalyse {
class CFieldStrengthAnalyse : public CSpreadAnalyse {
 public:
  IStations *pStations;
  float_t hm;
  string outputPath;
  IFileFloatArray *pData;
  float_t offsetDB;
  bool needComputeAll;
  vector<double_t> crsv;
  Station stationInfo;
  /// @brief
  /// 分析范围，默认是整个dem，每个子类有自己的get/putExtent方法，为该数组设置值
  double_t subExtent[4];

 public:
  CFieldStrengthAnalyse();
  virtual ~CFieldStrengthAnalyse() {
    if (pData != nullptr) {
      delete pData;
      pData = nullptr;
    }
    // if (pStations != nullptr) {
    //   delete pStations;
    //   pStations = nullptr;
    // }
  }
  // 场强分析算法的调用入口，每个场强分析类都有自己的场强分析调用入口，重定义。
  bool FieldStrengthAnalyse(string savePath, RasterCreateFileType type);

 protected:
  /* 子类必须实现的方法 */
  /// @brief 发射站覆盖分析
  /// @param stationInfo 台站信息
  /// @param rsv
  /// 动态数组，用于保存每个子类保存的数据，每次使用前需要留意数组中已经存在的数据
  /// @param point 待分析的点坐标
  /// @return 返回点对点之间场强分析的结果
  virtual float_t GetRadiuValue(const Station &stationInfo,
                                vector<double_t> *rsv,
                                const OGRPoint &point) = 0;

  /// @brief
  /// @param stationInfo 台站信息
  /// @param rsv
  /// @param point 待分析的点坐标
  /// @return
  virtual float_t GetRadiuValueRev(const Station &stationInfo,
                                   const vector<double_t> &rsv,
                                   const OGRPoint &point) const = 0;

  /// @brief
  /// @param stationInfo
  /// @param rsv
  /// 是一个数组，该数组里面装了台站的初始的发射功率以及其他相关信息
  virtual void PrepareReservedValues(const Station &stationInfo,
                                     vector<double_t> *rsv) = 0;

  /// @brief 分析所需要的其他栅格数据
  /// @return 返回准备的状态
  virtual bool PrepareOtherData() = 0;

  /// @brief 获取模型的名称
  /// @return 返回模型的名称
  virtual string GetModelName() = 0;

 protected:
  void ComputeOneStation(const Station &para);
  bool ComputeOneStationByExtent(const Station &para, vector<double_t> *rsv);
  int64_t *GetSubRowCol(double_t xmin, double_t ymin, double_t xmax,
                        double_t ymax);
};

}  // namespace fieldstrengthanalyse
}  // namespace spreadanalyse
}  // namespace spread

#endif  // INCLUDE_SPREAD_MODEL_FIELDSTRENGTHANALYSE_H_
