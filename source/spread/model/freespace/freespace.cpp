/**
 * @copyright ALL COPYRIGH RESERVED BY SHIGP
 * @author shigp
 */

#include <spread/model/freespace/freespace.h>

using spread::spreadbase::combine::freespace::CFreeSpace;

/* ---------------CFreeSpace自由空间传播模型的实现-------------------- */
void CFreeSpace::PrepareReservedValues(const Station& para,
                                       std::vector<double>* rsv) {
  rsv->clear();
  double dV = para.power + 30 + 95.23;
  rsv->push_back(dV);
}

/// @brief 自由传播算法的实现
/// @param stationInfo 站信息
/// @param rsv
/// @param point 要计算的点坐标
/// @return 返回point位置处的场强
float_t CFreeSpace::GetRadiuValue(const Station& stationInfo,
                                  std::vector<double_t>* rsv,
                                  const OGRPoint& point) {
  float_t tfV;
  double_t X = point.getX();
  double_t Y = point.getY();
  double_t Z = point.getZ();
  float_t d =
      sqrt(pow(stationInfo.x - X, 2.0f) + pow(stationInfo.y - Y, 2.0f) +
           pow(stationInfo.stationHeight + stationInfo.dem - Z - hm, 2.0f));
  tfV = 32.4 + 20 * log10(d / 1000) + 20 * log10(stationInfo.frequency);
  double_t dbLoss = CCombine::GetRadiuValue(stationInfo, rsv, point);
  // return rsv.ReservedValues[0] - tfV - dbLoss;
  return rsv->front() - tfV - dbLoss;
}

//
void CFreeSpace::FieldStrengthAnalyse(std::string outPutPath,
                                      RasterCreateFileType type, bool* pVal) {
  *pVal = CCombine::FieldStrengthAnalyse(outPutPath, type);
}
