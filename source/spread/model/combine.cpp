/**
 * @copyright ALL COPYRIGH RESERVED BY SHIGP
 * @author shigp
 */

#include <spread/model/combine.h>

#include <iostream>

using spatdata::PBand_T;
using spread::spreadbase::combine::CCombine;
using std::cout;
using std::endl;
using std::string;

bool CCombine::FieldStrengthAnalyse(std::string savePath,
                                    RasterCreateFileType type) {
  bool isOk = CFieldStrength::FieldStrengthAnalyse(savePath, type);
  return isOk;
}
// 重写接口方法
float_t CCombine::GetRadiuValue(const Station& stationInfo,
                                std::vector<double_t>* rsv,
                                const OGRPoint& point) {
  // 可能有多个绕射模型，所以els中可能有多个不同的绕射模型。
  double_t X = point.getX();
  double_t Y = point.getY();
  double_t Z = point.getZ();
  double_t allNum = 0;
  double_t loss = 0;
  for (int k = els.size() - 1; k >= 0; k--) {
    OGRPoint point;
    point.setX(X);
    point.setY(Y);
    point.setZ(Z);
    loss = els.at(k)->GetDBLoss(stationInfo, point);
    allNum += loss;
  }
  cout << rsv->size() << std::endl;
  return allNum;
}

string CCombine::GetModelName() { return "CombineAnalyse"; }

float_t CCombine::GetRadiuValueRev(const Station& para,
                                   const std::vector<double_t>& rsv,
                                   const OGRPoint& point) const {
  double AllNum = 0;
  for (int k = els.size() - 1; k >= 0; k--) {
    double loss;
    loss = els.at(k)->GetDBLossRev(para, point);
    AllNum += loss;
  }

  std::cout << rsv.size() << std::endl;
  return AllNum;
}
bool CCombine::PrepareOtherData() {
  for (int k = 0; k < static_cast<int>(otherPaths.size()); k++) {
    bool IsOk;
    PBand_T pBand;
    pBand.bandNumber = 1;
    OGRPoint* lt;
    otherReaders[k]->OpenRaster(otherPaths[k], &IsOk);
    if (!IsOk) {
      errorInfo = "打开" + otherNames[k] + "失败";
      return false;
    }
    otherReaders[k]->SetRasterBand(pBand, &IsOk);
    if (!IsOk) {
      errorInfo = "打开" + otherNames[k] + "失败";
      return false;
    }
    pEnvi->GetLeftTop(&lt);

    IFileFloatArray* pArray;
    otherReaders[k]->GetBlockDataByCoord(lt, cellSize, cols, rows, noData,
                                         &pArray);
    if (otherDatas.at(k) != nullptr) {
      std::vector<IFileFloatArray*>::iterator pos = otherDatas.begin() + k;
      otherDatas.erase(pos);
    }
    otherDatas.push_back(pArray);
    if (pArray == nullptr) {
      errorInfo = "打开" + otherNames[k] + "失败";
      return false;
    }
  }
  for (int k = els.size() - 1; k >= 0; k--) {
    els.at(k)->prepareAnalyseEnvi();
  }
  return true;
}
