/**
 * @copyright all Copyright reserved by shigp
 * @author shigp
 */

#ifndef INCLUDE_SPREAD_MODEL_FREESPACE_FREESPACE_H_
#define INCLUDE_SPREAD_MODEL_FREESPACE_FREESPACE_H_

#include <spread/model/combine.h>

#include <string>
#include <vector>

using spread::spreadbase::combine::CCombine;
using std::string;
using std::vector;

namespace spread {
namespace spreadbase {
namespace combine {
namespace freespace {

class CFreeSpace : public CCombine {
 protected:
  float_t GetRadiuValue(const Station &stationInfo, vector<double> *rsv,
                        const OGRPoint &point);
  void PrepareReservedValues(const Station &para, vector<double> *rsv);

 public:
  // 调用父类CFreeSpaceAnalyse的FieldStrengthAnalyse方法
  void FieldStrengthAnalyse(string savePath, RasterCreateFileType type,
                            bool *pVal);
};
}  // namespace freespace
}  // namespace combine
}  // namespace spreadbase
}  // namespace spread

#endif  // INCLUDE_SPREAD_MODEL_FREESPACE_FREESPACE_H_
