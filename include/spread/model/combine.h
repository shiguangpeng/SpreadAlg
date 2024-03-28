/**
 * @copyright ALL COPYRIGH RESERVED BY SHIGP
 * @author shigp
 */

#ifndef INCLUDE_SPREAD_MODEL_COMBINE_H_
#define INCLUDE_SPREAD_MODEL_COMBINE_H_

#include <spread/model/dbloss/dbloss.h>
#include <spread/model/fieldstrength.h>

#include <string>
#include <vector>

using spread::IDBLossElement;
using spread::spreadbase::fieldstrength::CFieldStrength;

namespace spread {
namespace spreadbase {
namespace combine {
/**
 * @brief Combine联合分析实现类
 */
class CCombine : public CFieldStrength {
 public:
  CCombine() = default;
  ~CCombine() = default;

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
}  // namespace spreadbase
}  // namespace spread

#endif  // INCLUDE_SPREAD_MODEL_COMBINE_H_
