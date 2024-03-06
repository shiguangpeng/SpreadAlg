#include <fmt/format.h>
#include <spread/spread.h>

#include <iostream>

#include "spread/spatidatamanager.h"

using namespace spread;
using namespace spatidatamanager;

/*** --------------------- 所有场强分析算法的基类，实现 ------------------------------- */

CSpreadAnalyse::CSpreadAnalyse() {}
CSpreadAnalyse::~CSpreadAnalyse() {}
bool CSpreadAnalyse::InitEnvironment() {
  // 通过打开的栅格数据，初始化类相关成员。
  const char* path = elevationPath.c_str();
  RasterDataSource ptr = RasterDataSource();
  GDALDataset* data = ptr.OpenRaster(path);
  if (data == nullptr) {
    errorInfo = "栅格数据打开失败。";
    return false;
  }
  // 使用GDAL数据指针初始化栅格属性
  // 行列数
  cols = data->GetRasterXSize();
  rows = data->GetRasterYSize();
  // 仿射变换参数
  double_t geoInfo[6];
  data->GetGeoTransform(geoInfo);
  xMin = geoInfo[0];
  yMax = geoInfo[3];
  cellSize = geoInfo[1];
  // 获取当前tif图的所有波段数量
  int bandNum = data->GetRasterCount();
  // 波段从1开始，dem这种结果图层也只有一个图层
  int* dataTypeFlag = nullptr;
  double_t bandNoData = 0.0;
  for (int i = 1; i <= bandNum; i++) {
    GDALRasterBand* band = data->GetRasterBand(i);
    GDALDataType dataType = band->GetRasterDataType();
    // 根据GetNoDataValue()方法的说明，特定数据类型需要使用其他方法
    if (dataType == GDALDataType::GDT_Int64) {
      bandNoData = band->GetNoDataValueAsInt64(dataTypeFlag);
    } else if (dataType == GDALDataType::GDT_UInt64) {
      bandNoData = band->GetNoDataValueAsUInt64(dataTypeFlag);
    } else {
      bandNoData = band->GetNoDataValue(dataTypeFlag);
    }
    // todo: 就项目来说，像dem这种结果tiff图，土地利用类型等，一般只有一个波段。这里先break;
    break;
  }
  if (dataTypeFlag) {
    noData = bandNoData;
  }
  // 结束前关闭数据集对象
  GDALClose(data);
  GDALDestroy();
  errorInfo = "栅格数据打开成功。";
  return true;
}

/** ---------------------CFieldStrengthAnalyse 场强分析类的实现 ---------------------------------*/
bool CFieldStrengthAnalyse::FieldStrengthAnalyse(std::string savePath, RasterCreateFileType type) {
  errorInfo = "";
  long Count;
  pStations->GetCount(&Count);
  if (Count == 0) {
    errorInfo = "参与计算的台站数为0";
    return false;
  }
  bool IsOk;
  if (!CSpreadAnalyse::InitEnvironment()) {
    return false;
  }

  if (!PrepareOtherData()) {
    return false;
  }

  // COM组件的写法，要将接口的引用计数减1，当减为0时，会使用delete将该接口（struct）释放
  // if (pData != NULL) {
  //   pData->Release();
  // }

  long Cols, Rows;
  pEnvi->GetCols(&Cols);
  pEnvi->GetRows(&Rows);
  pData->SetSize(Cols * Rows, noData, &IsOk);
  if (!IsOk) {
    errorInfo = "数组初始化失败";
    return false;
  }
  for (int k = 0; k < Count; k++) {
    Station para;
    pStations->GetItem(k, &para);
    ComputeOneStation(para);
  }
  // 保存结果，生成tiff
  // if (!SaveOutput(pData, type, RasterCreateDataType::rcdtFloat32)) {
  //   errorInfo = "保存数据失败, 可能是场强阈值太大覆盖范围为0或目标路径有误";
  //   return false;
  // }
  return true;
}

void CFieldStrengthAnalyse::ComputeOneStation(Station& para, std::string progress = "") {
  if (progress.empty()) {
    progress = para.name;
    progress = "计算" + GetModelName() + "场强-" + progress;
  }
  std::vector<double_t> rsv;
  PrepareReservedValues(para, rsv);
  ComputeOneStationByExtent(para, rsv);
}

bool CFieldStrengthAnalyse::ComputeOneStationByExtent(Station& para, std::vector<double_t>& rsv) {
  //
  int MaxRadius = std::max(cols, rows);
  // 计算当前分析点的行列号
  int Col, Row;
  Col = (para.x - xMin) / cellSize;
  if (Col < 0)
    Col = 0;
  else if (Col >= cols)
    Col = cols - 1;
  Row = (yMax - para.y) / cellSize;
  if (Row < 0)
    Row = 0;
  else if (Row >= rows)
    Row = rows - 1;
  // 分析范围默认为DEM范围
  int tempextent[4] = {0, cols, 0, rows};
  int* pextent = tempextent;
  int sr = 0;
  int er = MaxRadius;
  // 计算分析范围所属行列号范围
  if (subExtent[0] != 0 || subExtent[1] != 0 || subExtent[2] != 0 || subExtent[3] != 0) {
    pextent = GetSubRowCol(subExtent[0], subExtent[1], subExtent[2], subExtent[3]);
    // 计算分析点到分析范围内移动的最小与最大行列数
    int* movecr = new int[4];
    movecr[0] = abs(Col - pextent[0]);
    movecr[1] = abs(Col - pextent[1]);
    movecr[2] = abs(Row - pextent[2]);
    movecr[3] = abs(Row - pextent[3]);
    // 计算移动的起止行列索引
    sr = movecr[0];
    er = movecr[0];
    for (int i = 0; i < 4; i++) {
      if (movecr[i] < sr) {
        sr = movecr[i];
      }
      if (movecr[i] > er) {
        er = movecr[i];
      }
    }
    // 如果台站位于分析范围内则起止行列索引为0
    if (para.x - subExtent[0] >= 0 && para.y - subExtent[1] >= 0 && para.x - subExtent[2] <= 0
        && para.y - subExtent[3] <= 0) {
      sr = 0;
    }
  }
  CRect ComputeRect(pextent[0], pextent[2], pextent[1], pextent[3]);

  if (needComputeAll) {
    sr = 0;
    er = MaxRadius;
    ComputeRect.left = 0;
    ComputeRect.right = cols;
    ComputeRect.top = 0;
    ComputeRect.bottom = rows;
    pextent[0] = 0;
    pextent[1] = cols;
    pextent[2] = 0;
    pextent[3] = rows;
  }
  // 从最小行列号开始分析
  for (int Radius = sr; Radius < er; Radius++) {
    bool HasCompute = false;
    if ((Col - Radius < 0) && (Col + Radius >= cols) && (Row - Radius < 0)
        && (Row + Radius >= rows)) {
      break;
    }
    double_t X, Y;
    float_t Z;
    int i, j;
    long Posi;
    float_t result;

    long nodata = noData;
    // 上
    i = Row - Radius;
    if (i >= 0 && i >= pextent[2] && i <= pextent[3]) {
      Y = yMax - cellSize * i - cellSize / 2;
      int FromX = ComputeRect.left;
      if (FromX < 0) FromX = 0;
      int ToX = ComputeRect.right;
      if (ToX >= cols) ToX = cols - 1;
      X = xMin + FromX * cellSize + cellSize / 2;
      Posi = FromX + i * cols;
      for (j = FromX; j <= ToX; j++) {
        pElevData->GetValueAsFloat(Posi, &Z);
        if (nodata == (long)Z) {
          X += cellSize;
          Posi++;
          continue;
        }
        OGRPoint p(X, Y, Z);
        result = GetRadiuValue(para, rsv, p) + offsetDB;
        if (nodata == (long)result) {
          X += cellSize;
          Posi++;
          continue;
        }
        if (result >= para.freqThreshold) {
          float FormerValue;
          pData->GetValueAsFloat(Posi, &FormerValue);
          if ((long)FormerValue == (long)noData) {
            pData->SetValueAsFloat(Posi, result);
          } else if (result > FormerValue) {
            pData->SetValueAsFloat(Posi, result);
          }
          HasCompute = true;
        }
        X += cellSize;
        Posi++;
      }
    }
    // 下
    i = Row + Radius;
    if (i < rows && i <= pextent[3] && i >= pextent[2]) {
      Y = yMax - cellSize * i - cellSize / 2;
      int FromX = ComputeRect.left;
      if (FromX < 0) FromX = 0;
      int ToX = ComputeRect.right;
      if (ToX >= cols) ToX = cols - 1;
      X = xMin + FromX * cellSize + cellSize / 2;
      Posi = FromX + i * cols;
      for (j = FromX; j <= ToX; j++) {
        pElevData->GetValueAsFloat(Posi, &Z);
        if (nodata == (long)Z) {
          X += cellSize;
          Posi++;
          continue;
        }
        OGRPoint point(X, Y, Z);
        result = GetRadiuValue(para, rsv, point) + offsetDB;
        if (nodata == (long)result) {
          X += cellSize;
          Posi++;
          continue;
        }
        if (result >= para.freqThreshold) {
          float FormerValue;
          pData->GetValueAsFloat(Posi, &FormerValue);
          if (nodata == (long)FormerValue) {
            pData->SetValueAsFloat(Posi, result);
          } else if (result > FormerValue) {
            pData->SetValueAsFloat(Posi, result);
          }
          HasCompute = true;
        }
        X += cellSize;
        Posi++;
      }
    }
    // 左
    j = Col - Radius;
    if (j >= 0 && j >= pextent[0] && j <= pextent[1]) {
      X = xMin + cellSize * j + cellSize / 2;
      int FromY = ComputeRect.top + 1;
      if (FromY < 0) FromY = 0;
      int ToY = ComputeRect.bottom - 1;
      if (ToY >= rows) ToY = rows - 1;
      Y = yMax - FromY * cellSize - cellSize / 2;
      Posi = FromY * cols + j;
      for (i = FromY; i <= ToY; i++) {
        pElevData->GetValueAsFloat(Posi, &Z);
        if (nodata == (long)Z) {
          Y -= cellSize;
          Posi += cols;
          continue;
        }
        OGRPoint point(X, Y, Z);
        result = GetRadiuValue(para, rsv, point) + offsetDB;
        if ((long)result == (long)noData) {
          Y -= cellSize;
          Posi += cols;
          continue;
        }
        if (result >= para.freqThreshold) {
          float FormerValue;
          pData->GetValueAsFloat(Posi, &FormerValue);
          if (nodata == (long)FormerValue) {
            pData->SetValueAsFloat(Posi, result);
          } else if (result > FormerValue) {
            pData->SetValueAsFloat(Posi, result);
          }
          HasCompute = true;
        }
        Y -= cellSize;
        Posi += cols;
      }
    }
    // 右
    j = Col + Radius;
    if (j < cols && j <= pextent[1] && j >= pextent[0]) {
      X = xMin + cellSize * j + cellSize / 2;
      int FromY = ComputeRect.top + 1;
      if (FromY < 0) FromY = 0;
      int ToY = ComputeRect.bottom - 1;
      if (ToY >= rows) ToY = rows - 1;
      Y = yMax - FromY * cellSize - cellSize / 2;
      Posi = FromY * cols + j;
      for (i = FromY; i <= ToY; i++) {
        pElevData->GetValueAsFloat(Posi, &Z);
        if (nodata == (long)Z) {
          Y -= cellSize;
          Posi += cols;
          continue;
        }
        OGRPoint point(X, Y, Z);
        result = GetRadiuValue(para, rsv, point) + offsetDB;
        if ((long)result == (long)noData) {
          Y -= cellSize;
          Posi += cols;
          continue;
        }
        if (result >= para.freqThreshold) {
          float FormerValue;
          pData->GetValueAsFloat(Posi, &FormerValue);
          if (nodata == (long)FormerValue) {
            pData->SetValueAsFloat(Posi, result);
          } else if (result > FormerValue) {
            pData->SetValueAsFloat(Posi, result);
          }
          HasCompute = true;
        }
        Y -= cellSize;
        Posi += cols;
      }
    }

    if (!needComputeAll && !HasCompute) {
      break;
    }
  }
  return true;
}

int* CFieldStrengthAnalyse::GetSubRowCol(double xmin, double ymin, double xmax, double ymax) {
  // 注意：在堆上分配了内存，注意释放
  int* result = new int[4];
  int colmin, rowmin, colmax, rowmax;
  colmin = (xmin - xMin) / cellSize;
  colmax = (xmax - xMin) / cellSize;
  rowmin = (yMax - ymax) / cellSize;
  rowmax = (yMax - ymin) / cellSize;

  result[0] = colmin;
  result[1] = colmax;
  result[2] = rowmin;
  result[3] = rowmax;
  return result;
}

// bool CFieldStrengthAnalyse::SaveOutput(IFileFloatArray* pData, RasterCreateFileType type,
//                                        RasterCreateDataType dtype) {
//   int l = cols - 1, t = rows - 1, r = 0, b = 0;
//   long pos = 0;
//   long nodata = (long)noData;
//   for (int i = 0; i < rows; i++) {
//     for (int j = 0; j < cols; j++) {
//       float fV;
//       pData->GetValueAsFloat(pos, &fV);
//       if (nodata != (long)fV) {
//         if (l > j) l = j;
//         if (r < j) r = j;
//         if (t > i) t = i;
//         if (b < i) b = i;
//       }
//       pos++;
//     }
//   }
//   if ((l > r) || (t > b)) return false;
//   IPoint* LeftTop;
//   ::CoCreateInstance(CLSID_PointDef, NULL, CLSCTX_INPROC_SERVER, IID_IPoint, (LPVOID*)&LeftTop);
//   LeftTop->put_X(XMin + l * CellSize);
//   LeftTop->put_Y(YMax - t * CellSize);
//   int rCols;
//   int rRows;
//   double_t OutputCellSize;
//   pEnvi->get_OutputCellSize(&OutputCellSize);
//   rCols = (r * CellSize - l * CellSize + CellSize) / OutputCellSize;
//   rRows = (b * CellSize - t * CellSize + CellSize) / OutputCellSize;
//   IGDALRasterCreator* pCreator;
//   ISpatialReference* pSpatial;
//   IGDALRasterProperties* pPro;
//   pElevs->QueryInterface(IID_IGDALRasterProperties, (void**)&pPro);
//   pPro->Release();
//   pPro->get_SpatialReference(&pSpatial);
//   ::CoCreateInstance(CLSID_GDALRasterCreator, NULL, CLSCTX_INPROC_SERVER, IID_IGDALRasterCreator,
//                      (LPVOID*)&pCreator);
//   if (progress != NULL) progress->BeginProgress(CComBSTR("创建结果图层"));
//   VARIANT_BOOL IsOk;
//   pCreator->CreateRaster(CComBSTR(OutPath), LeftTop, rCols, rRows, OutputCellSize, type, dtype,
//   1,
//                          pSpatial, NoData, &IsOk);
//   pCreator->Release();
//   pSpatial->Release();
//   LeftTop->Release();
//   if (!IsOk) {
//     ErrorInfo = "创建数据失败";
//     return false;
//   }
//   IGDALRasterSaverByFileArray* pSaver;
//   ::CoCreateInstance(CLSID_GDALRasterSaverByFileArray, NULL, CLSCTX_INPROC_SERVER,
//                      IID_IGDALRasterSaverByFileArray, (LPVOID*)&pSaver);
//   pSaver->OpenRaster(CComBSTR(OutPath), &IsOk);
//   if (!IsOk) {
//     ErrorInfo = "保存数据失败";
//     pSaver->Release();
//     return false;
//   }
//   VARIANT pBand;
//   pBand.vt = VT_I4;
//   pBand.intVal = 1;
//   pSaver->PutRasterBand(pBand, &IsOk);
//   if (!IsOk) {
//     ErrorInfo = "保存数据失败";
//     pSaver->Release();
//     return false;
//   }
//   if ((dtype == rcdtFloat32) || (dtype == rcdtFloat64))
//     pSaver->SaveBlockDataWithInterpolateByRect(0, 0, rCols - 1, rRows - 1, l, t, r, b, Cols,
//     Rows,
//                                                pData, progress, &IsOk);
//   else
//     pSaver->SaveBlockDataByRect(0, 0, rCols - 1, rRows - 1, l, t, r, b, Cols, Rows, pData,
//     progress,
//                                 &IsOk);
//   if (!IsOk) {
//     ErrorInfo = "保存数据失败";
//     pSaver->Release();
//     return false;
//   }
//   pSaver->Release();
//   return true;
// }

/*** ----------------------CCombineAnalyse 混合传播模型的实现------------------------------ */

float_t CCombineAnalyse::GetRadiuValue(Station& stationInfo, std::vector<double_t>& rsv,
                                       OGRPoint* point) {
  // 可能有多个绕射模型，所以els中可能有多个不同的绕射模型。
  double_t X = point->getX();
  double_t Y = point->getY();
  double_t Z = point->getZ();
  double_t allNum = 0;
  for (int k = els.size() - 1; k >= 0; k--) {
    double_t loss;
    els.at(k)->GetDBLoss(stationInfo, X, Y, Z, &loss);
    allNum += loss;
  }
  return allNum;
}

bool CCombineAnalyse::FieldStrengthAnalyse(std::string savePath, RasterCreateFileType type) {
  bool isOk = CFieldStrengthAnalyse::FieldStrengthAnalyse(savePath, type);
  return isOk;
}

/*** ----------------------CFreeSpaceAnalyse 自由空间传播模型的实现------------------------------ */

/// @brief 自由传播算法的实现
/// @param stationInfo 站信息
/// @param rsv
/// @param point 要计算的点坐标
/// @return 返回point位置处的场强
float_t CFreeSpaceAnalyse::GetRadiuValue(Station& stationInfo, std::vector<double_t>& rsv,
                                         OGRPoint* point) {
  float_t tfV;
  double_t X = point->getX();
  double_t Y = point->getY();
  double_t Z = point->getZ();
  float d = sqrt(pow(stationInfo.x - X, (double)2.0) + pow(stationInfo.y - Y, (double)2.0)
                 + pow(stationInfo.stationHeight + stationInfo.dem - Z - hm, (double)2.0));
  tfV = 32.4 + 20 * log10(d / 1000) + 20 * log10(stationInfo.frequency);
  OGRPoint point(X, Y, Z);
  double_t dbLoss = CCombineAnalyse::GetRadiuValue(stationInfo, rsv, point);
  // return rsv.ReservedValues[0] - tfV - dbLoss;
  return rsv.front() - tfV - dbLoss;
}