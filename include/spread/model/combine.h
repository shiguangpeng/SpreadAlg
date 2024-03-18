/**
 * @copyright ALL COPYRIGH RESERVED BY SHIGP
 * @author shigp
 */

#ifndef INCLUDE_SPREAD_MODEL_COMBINE_H_
#define INCLUDE_SPREAD_MODEL_COMBINE_H_

#include <spread/model/dbloss/dbloss.h>
#include <spread/model/fieldstrengthanalyse.h>

#include <string>
#include <vector>

using dbloss::IDBLossElement;
using spread::spreadanalyse::fieldstrengthanalyse::CFieldStrengthAnalyse;

namespace spread {
namespace spreadanalyse {
namespace combine {
/**
 * @brief 联合分析类的接口
 */
// struct ICombineAnalyse {
//   virtual void AddOtherRaster(string path);
//   virtual void GetOtherRaster(int64_t nIndex, IFileFloatArray **pVal);
//   virtual void get_ElevationData(IFileFloatArray **pVal);
//   virtual void get_Elevation(IGDALRasterReaderByFileArray **pVal);
//   virtual void get_NoData(double *pVal);
//   virtual void AddDBLossElement(IDBLossElement *el);
// };

class CCombineAnalyse : public CFieldStrengthAnalyse {
 public:
  CCombineAnalyse() = default;
  ~CCombineAnalyse() = default;

 public:
  bool FieldStrengthAnalyse(string savePath, RasterCreateFileType type);

 protected:
  float_t GetRadiuValue(const Station &stationInfo, vector<double_t> *rsv,
                        const OGRPoint &point) override;

  float_t GetRadiuValueRev(const Station &para, const vector<double_t> &rsv,
                           const OGRPoint &point) const override;

  bool PrepareOtherData() override;
  std::string GetModelName() override;

 protected:
  vector<IGDALRasterReaderByFileArray *> otherReaders;
  vector<IFileFloatArray *> otherDatas;
  vector<IDBLossElement *> els;
  vector<string> otherNames;
  vector<string> otherPaths;
};

}  // namespace combine
}  // namespace spreadanalyse
}  // namespace spread

#endif  // INCLUDE_SPREAD_MODEL_COMBINE_H_
