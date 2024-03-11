#include <fmt/format.h>
#include <spread/spread.h>

#include <iostream>
#include <vector>

using namespace spread;

/* -----------------------------各个类的析构函数实现------------------------------------*/
// CSpreadAnalyse::~CSpreadAnalyse() = default;
CSpreadAnalyse::CSpreadAnalyse() {
  pEnvi = nullptr;
  errorInfo = "";
  // 初始化pElevs对象
  pElevs = nullptr;
  pElevData = nullptr;
  isInit = false;
  otherDataPreload = false;
}

/*** --------------------- 所有场强分析算法的基类，实现 ------------------------------- */
bool CSpreadAnalyse::InitEnvironment() {
  // // 通过打开的栅格数据，初始化类相关成员。
  // const char* path = elevationPath.c_str();
  // RasterDataSource ptr = RasterDataSource();
  // GDALDataset* data = ptr.OpenRaster(path);
  // if (data == nullptr) {
  //   errorInfo = "栅格数据打开失败。";
  //   return false;
  // }
  // // 使用GDAL数据指针初始化栅格属性
  // // 行列数
  // cols = data->GetRasterXSize();
  // rows = data->GetRasterYSize();
  // // 仿射变换参数
  // double_t geoInfo[6];
  // data->GetGeoTransform(geoInfo);
  // xMin = geoInfo[0];
  // yMax = geoInfo[3];
  // cellSize = geoInfo[1];
  // // 获取当前tif图的所有波段数量
  // int bandNum = data->GetRasterCount();
  // // 波段从1开始，dem这种结果图层也只有一个图层
  // int* dataTypeFlag = nullptr;
  // double_t bandNoData = 0.0;
  // for (int i = 1; i <= bandNum; i++) {
  //   GDALRasterBand* band = data->GetRasterBand(i);
  //   GDALDataType dataType = band->GetRasterDataType();
  //   // 根据GetNoDataValue()方法的说明，特定数据类型需要使用其他方法
  //   if (dataType == GDALDataType::GDT_Int64) {
  //     bandNoData = band->GetNoDataValueAsInt64(dataTypeFlag);
  //   } else if (dataType == GDALDataType::GDT_UInt64) {
  //     bandNoData = band->GetNoDataValueAsUInt64(dataTypeFlag);
  //   } else {
  //     bandNoData = band->GetNoDataValue(dataTypeFlag);
  //   }
  //   // todo: 就项目来说，像dem这种结果tiff图，土地利用类型等，一般只有一个波段。这里先break;
  //   break;
  // }
  // if (dataTypeFlag) {
  //   noData = bandNoData;
  // }
  // // 结束前关闭数据集对象
  // GDALClose(data);
  // GDALDestroy();
  // errorInfo = "栅格数据打开成功。";
  // return true;
  if (isInit) return true;
  bool IsOk;
  pElevs->OpenRaster(elevationPath, &IsOk);
  if (!IsOk) {
    errorInfo = "打开高程数据失败";
    return false;
  }
  // VARIANT pBand;
  // pBand.vt = VT_I4;
  // pBand.intVal = 1;
  int pBand = 1;
  // 设置栅格的波段，默认第一个波段
  PBand_T band;
  band.bandNumber = pBand;
  pElevs->SetRasterBand(band, &IsOk);
  if (!IsOk) {
    errorInfo = "打开高程数据失败";
    return false;
  }
  IGDALRasterProperties* pPro = nullptr;
  // pElevs->QueryInterface(IID_IGDALRasterProperties, (void**)&pPro);
  pPro->GetNoData(&noData);
  pPro->GetCols(&cols);
  pPro->GetRows(&rows);
  OGREnvelope* Extent;
  pPro->GetExtent(&Extent);
  OGRPoint* ppt = nullptr;
  // ::CoCreateInstance(CLSID_PointDef, nullptr, CLSCTX_INPROC_SERVER, IID_IPoint, (void**)&ppt);
  double_t left, top;
  // Extent->GetLeft(&left);
  // Extent->GetTop(&top);
  left = Extent->MinY;
  top = Extent->MaxX;
  // ppt->PutCoord(left, top);
  ppt->setX(top);
  ppt->setY(left);
  pPro->GetCellSize(&cellSize);
  if (pEnvi == nullptr) {
    // ::CoCreateInstance(CLSID_AnalyseEnvi, nullptr, CLSCTX_INPROC_SERVER, IID_IAnalyseEnvi,
    //                    (LPVOID*)&pEnvi);  // 得到接口指针
    pEnvi->SetLeftTop(ppt);
    pEnvi->SetCols(cols);
    pEnvi->SetRows(rows);
    pEnvi->SetCellSize(cellSize);
    pEnvi->SetOutputCellSize(cellSize);
  } else {
    double_t rCellSize;
    pEnvi->GetOutputCellSize(&rCellSize);
    if (rCellSize == 0) {
      pEnvi->SetOutputCellSize(cellSize);
    }
    pEnvi->GetCellSize(&rCellSize);
    if (rCellSize == 0)
      pEnvi->SetCellSize(cellSize);
    else
      cellSize = rCellSize;
    OGRPoint* lt;
    pEnvi->GetLeftTop(&lt);
    if (lt != nullptr) {
      ppt = lt;
    } else
      pEnvi->SetLeftTop(ppt);
    long rRows;
    long rCols;
    pEnvi->GetRows(&rRows);
    pEnvi->GetCols(&rCols);
    if (rCols == 0)
      pEnvi->SetCols(cols);
    else
      cols = rCols;
    if (rRows == 0)
      pEnvi->SetRows(rows);
    else
      rows = rRows;
  }
  xMin = ppt->getX();
  yMax = ppt->getY();
  pElevData = nullptr;
  OGRPoint* lt;
  pEnvi->GetLeftTop(&lt);
  pElevs->GetBlockDataByCoord(lt, cellSize, cols, rows, noData, &pElevData);
  isInit = true;
  return true;
}
/*-----------------------------------CRect的实现-----------------------------------------------*/
CRect::CRect() = default;
CRect::CRect(int left, int top, int right, int bottom) {
  this->left = left;
  this->top = top;
  this->bottom = bottom;
  this->right = right;
}

/** ---------------------CFieldStrengthAnalyse 场强分析类的实现 ---------------------------------*/
CFieldStrengthAnalyse::CFieldStrengthAnalyse() {
  hm = 2.5;
  pData = nullptr;
  offsetDB = 0;
  needComputeAll = false;
  subExtent[0] = 0;
  subExtent[1] = 0;
  subExtent[2] = 0;
  subExtent[3] = 0;
}
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
  // if (pData != nullptr) {
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
  // 为了避免gcc未使用变量检测
  std::cout << "inputpath: " << savePath << std::endl;
  std::cout << "output raster type: " << type << std::endl;
  return true;
}

void CFieldStrengthAnalyse::ComputeOneStation(Station& para) {
  std::cout << "计算" + GetModelName() + "场强-" << std::endl;

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
  long int tempextent[4] = {0, cols, 0, rows};
  long int* pextent = tempextent;
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

long int* CFieldStrengthAnalyse::GetSubRowCol(double xmin, double ymin, double xmax, double ymax) {
  // 注意：在堆上分配了内存，注意释放
  long int* result = new long int[4];
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
//   ::CoCreateInstance(CLSID_PointDef, nullptr, CLSCTX_INPROC_SERVER, IID_IPoint,
//   (LPVOID*)&LeftTop); LeftTop->put_X(XMin + l * CellSize); LeftTop->put_Y(YMax - t * CellSize);
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
//   ::CoCreateInstance(CLSID_GDALRasterCreator, nullptr, CLSCTX_INPROC_SERVER,
//   IID_IGDALRasterCreator,
//                      (LPVOID*)&pCreator);
//   if (progress != nullptr) progress->BeginProgress(CComBSTR("创建结果图层"));
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
//   ::CoCreateInstance(CLSID_GDALRasterSaverByFileArray, nullptr, CLSCTX_INPROC_SERVER,
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

bool CCombineAnalyse::FieldStrengthAnalyse(std::string savePath, RasterCreateFileType type) {
  bool isOk = CFieldStrengthAnalyse::FieldStrengthAnalyse(savePath, type);
  return isOk;
}
// 重写接口方法
float_t CCombineAnalyse::GetRadiuValue(Station& stationInfo, std::vector<double_t>& rsv,
                                       OGRPoint& point) {
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
  std::cout << rsv.size() << std::endl;
  return allNum;
}

std::string CCombineAnalyse::GetModelName() { return "CombineAnalyse"; }

float_t CCombineAnalyse::GetRadiuValueRev(Station& para, std::vector<double_t>& rsv,
                                          OGRPoint& point) {
  double AllNum = 0;
  for (int k = els.size() - 1; k >= 0; k--) {
    double loss;
    loss = els.at(k)->GetDBLossRev(para, point);
    AllNum += loss;
  }

  std::cout << rsv.size() << std::endl;
  return AllNum;
}
bool CCombineAnalyse::PrepareOtherData() {
  for (int k = 0; k < (int)otherPaths.size(); k++) {
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
    otherReaders[k]->GetBlockDataByCoord(lt, cellSize, cols, rows, noData, &pArray);
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

/*** ----------------------CFreeSpaceAnalyse 自由空间传播模型的实现------------------------------ */
void CFreeSpaceAnalyse::PrepareReservedValues(Station& para, std::vector<double>& rsv) {
  rsv.clear();
  double dV = para.power + 30 + 95.23;
  rsv.push_back(dV);
}

/// @brief 自由传播算法的实现
/// @param stationInfo 站信息
/// @param rsv
/// @param point 要计算的点坐标
/// @return 返回point位置处的场强
float_t CFreeSpaceAnalyse::GetRadiuValue(Station& stationInfo, std::vector<double_t>& rsv,
                                         OGRPoint& point) {
  float_t tfV;
  double_t X = point.getX();
  double_t Y = point.getY();
  double_t Z = point.getZ();
  float d = sqrt(pow(stationInfo.x - X, (double)2.0) + pow(stationInfo.y - Y, (double)2.0)
                 + pow(stationInfo.stationHeight + stationInfo.dem - Z - hm, (double)2.0));
  tfV = 32.4 + 20 * log10(d / 1000) + 20 * log10(stationInfo.frequency);
  double_t dbLoss = CCombineAnalyse::GetRadiuValue(stationInfo, rsv, point);
  // return rsv.ReservedValues[0] - tfV - dbLoss;
  return rsv.front() - tfV - dbLoss;
}

//
void CFreeSpaceAnalyse::FieldStrengthAnalyse(std::string outPutPath, RasterCreateFileType type,
                                             bool* pVal) {
  *pVal = CCombineAnalyse::FieldStrengthAnalyse(outPutPath, type);
}

// void CFreeSpaceAnalyse::FieldStrengthAnalyse(std::string savePath, RasterCreateFileType type,
//                                              bool* pVal) {
//   *pVal = CCombineAnalyse::FieldStrengthAnalyse(savePath, type);
// }

/*--------------------------AnalyseEnvironment---------------------------------------*/
AnalyseEnvironment::AnalyseEnvironment() {
  // 初始化类的成员变量
  leftTop = nullptr;
  cellSize = 0.0;
  cols = rows = 0;
  outPutCellSize = 0.0;
}
// 重写接口对成员属性的get/set方法，对成员变量获取或设置值
void AnalyseEnvironment::GetLeftTop(OGRPoint** point) { *point = leftTop; }
void AnalyseEnvironment::SetLeftTop(OGRPoint* point) {
  leftTop->setX(point->getX());
  leftTop->setY(point->getY());
  leftTop->setZ(point->getZ());
}

void AnalyseEnvironment::GetCellSize(double_t* cellSize) { *cellSize = this->cellSize; }
void AnalyseEnvironment::SetCellSize(double_t cellSize) { this->cellSize = cellSize; }

void AnalyseEnvironment::GetCols(long* cols) { *cols = this->cols; }
void AnalyseEnvironment::SetCols(long cols) { this->cols = cols; }

void AnalyseEnvironment::GetRows(long* rows) { *rows = this->rows; }
void AnalyseEnvironment::SetRows(long rows) { this->rows = rows; }

void AnalyseEnvironment::GetOutputCellSize(double_t* outPutCellSize) {
  *outPutCellSize = this->outPutCellSize;
}
void AnalyseEnvironment::SetOutputCellSize(double_t outPutCellSize) {
  this->outPutCellSize = outPutCellSize;
}

/* -------------------------------CStations的类定义-----------------------------------*/
void CStations::AddStation(Station station) { this->stations.push_back(&station); }
void CStations::GetCount(long* pVal) { *pVal = this->stations.size(); }
void CStations::RemoveAt(long index) {
  // 获取数组开始的迭代器
  std::vector<spread::Station*>::iterator spec = this->stations.begin() + index;
  this->stations.erase(spec);
}

void CStations::RemoveAll(void) { this->stations.clear(); }

void CStations::GetItem(long index, Station* pVal) {
  // std::vector<spread::Station*>::iterator spec = this->stations.begin() + index;
  *pVal = *this->stations.at(index);
}

void CStations::SetItem(long index, Station newVal) {
  std::vector<spread::Station*>::iterator spec = this->stations.begin() + index;
  this->stations.insert(spec, &newVal);
}

void CStations::PickHeightFromDEM(IGDALRasterReaderByPixel* pReader, bool* pVal) {
  int Size = stations.size();
  IGDALRasterProperties* pPro = nullptr;
  double CellSize;
  pPro->GetCellSize(&CellSize);
  OGREnvelope* rect;
  pPro->GetExtent(&rect);
  //
  // double XMin, YMin, XMax, YMax;
  double XMin, YMax;
  XMin = rect->MinX;
  // YMin = rect->MaxY;
  // XMax = rect->MaxX;
  YMax = rect->MaxY;
  long rows, cols;
  pPro->GetRows(&rows);
  pPro->GetCols(&cols);

  // 获取指定栅格的dem值
  float_t dV;
  for (int k = 0; k < Size; k++) {
    Station* para = stations.at(k);
    int Col = (para->x - XMin) / CellSize;
    int Row = (YMax - para->y) / CellSize;
    if ((Col < 0) || (Col >= cols) || (Row < 0) || (Row >= rows))
      para->dem = 0;
    else {
      pReader->GetPixelValue(Col, Row, &dV);
      para->dem = dV;
    }
  }
  // else {
  //   for (int k = 0; k < Size; k++) {
  //     Station* para = stations.at(k);
  //     int Col = (para->x - XMin) / CellSize;
  //     int Row = (YMax - para->y) / CellSize;
  //     if ((Col < 0) || (Col >= cols) || (Row < 0) || (Row >= rows))
  //       para->dem = 0;
  //     else {
  //       float dem;
  //       pReader->GetPixelValue(Col, Row, &dem);
  //       para->dem = dem;
  //     }
  //     pProgress->SetPos((float)k / Size * 100);
  //   }
  // }
  *pVal = true;
}

/** ----------------------CGDALRasterReaderByPixel的定义------------------------*/
GDALColorTable* CGDALRasterReaderByPixel::CreateColorTable() {
  GDALColorTable* pColorTable = nullptr;
  if (poBand == nullptr) {
    return nullptr;
  }
  GDALColorTable* gct = poBand->GetColorTable();
  if (gct == nullptr) {
    return nullptr;
  }
  pColorTable = new GDALColorTable(GPI_RGB);

  int colorKinds = gct->GetColorEntryCount();
  for (int k = 0; k < colorKinds; k++) {
    GDALColorEntry gce;
    gct->GetColorEntryAsRGB(k, &gce);
    pColorTable->SetColorEntry(k, &gce);
    // RGBColor nc;
    // nc.r = gce.c1;
    // nc.g = gce.c2;
    // nc.b = gce.c3;
    // pColorTable->AddColor(nc);
  }
  return pColorTable;
}
GDALRasterBand* CGDALRasterReaderByPixel::GetGDALBand() { return this->poBand; }
void CGDALRasterReaderByPixel::OpenRaster(std::string lpszPathName, bool* pVal) {
  if (poDataset != nullptr) {
    delete poDataset;
  }
  if (poSubDataset != nullptr) {
    delete poSubDataset;
  }
  poDataset = nullptr;
  poSubDataset = nullptr;
  rows = 0;
  cols = 0;
  cellSize = 0;
  poBand = nullptr;
  // extent->PutCoord(0, 0, 0, 0);
  extent->MaxX = 0;
  extent->MaxY = 0;
  extent->MinX = 0;
  extent->MinY = 0;
  metadata.clear();

  // 打开栅格数据文件，获得Dataset对象
  poDataset = (GDALDataset*)GDALOpen(lpszPathName.c_str(), GA_ReadOnly);
  if (poDataset == nullptr) {
    *pVal = false;
    return;
  }
  // const char* papszMetadata = GDALGetDriverShortName((GDALDriverH)poDataset);
  // The SUBDATASETS domain holds a list of child datasets. Normally this is used to provide
  // pointers to a list of images stored within a single multi image file.
  char** SUBDATASETS = GDALGetMetadata((GDALDatasetH)poDataset, "SUBDATASETS");
  if (CSLCount(SUBDATASETS) == 0) {
    rows = poDataset->GetRasterYSize();
    cols = poDataset->GetRasterXSize();
    double adfGeoTransform[6];
    poDataset->GetGeoTransform(adfGeoTransform);
    if (adfGeoTransform[5] > 0) {
      adfGeoTransform[3] = adfGeoTransform[3] + adfGeoTransform[5] * rows;
      adfGeoTransform[5] = -adfGeoTransform[5];
    }
    extent->MinX = adfGeoTransform[0];
    extent->MaxX = adfGeoTransform[0] + cols * adfGeoTransform[1];
    extent->MinY = adfGeoTransform[3] + rows * adfGeoTransform[5];
    extent->MaxY = adfGeoTransform[3];
    // extent->PutCoord(adfGeoTransform[0], adfGeoTransform[3],
    //                  adfGeoTransform[0] + cols * adfGeoTransform[1],
    //                  adfGeoTransform[3] + adfGeoTransform[5] * rows);
    cellSize = fabs(adfGeoTransform[1]);
  } else {
    for (int i = 0; SUBDATASETS[i] != nullptr; i++) {
      std::string tmpstr = SUBDATASETS[i];
      // 将每个SUBDATASETS按照等号分割
      int position = tmpstr.find("=", 0);
      // tmpstr = tmpstr.Right(tmpstr.GetLength() - tmpstr.Find("=") - 1);
      // 取等号右边的属性值
      metadata.push_back(tmpstr.substr(position));
    }
  }
  *pVal = true;
}
void CGDALRasterReaderByPixel::GetPixelValue(long col, long row, float_t* data) {
  if (poBand == nullptr) {
    *data = -INT_MAX;
    return;
  }
  CPLErr errorNumber = poBand->RasterIO(GF_Read, col, row, 1, 1, data, 1, 1, GDT_Float32, 0, 0);
  std::cout << "获取结果标志：" << errorNumber << std::endl;
}
void CGDALRasterReaderByPixel::SetRasterBand(PBand_T pBand, bool* pVal) {
  poBand = nullptr;
  if (poDataset == nullptr) {
    *pVal = false;
    return;
  }
  if (metadata.size() == 0) {
    poBand = poDataset->GetRasterBand(pBand.bandNumber);
  } else {
    rows = 0;
    cols = 0;
    cellSize = 0;
    // extent->PutCoord(0, 0, 0, 0);
    extent->MaxX = 0;
    extent->MaxY = 0;
    extent->MinX = 0;
    extent->MinY = 0;

    if (poSubDataset != nullptr) {
      delete poSubDataset;
    }
    poSubDataset = nullptr;
    // _bstr_t path = vr.bstrVal;
    // char* cpath = path;

    poSubDataset = (GDALDataset*)GDALOpen(pBand.rasterPath, GA_ReadOnly);
    if (poSubDataset == nullptr) {
      *pVal = false;
      return;
    }
    rows = poSubDataset->GetRasterYSize();
    cols = poSubDataset->GetRasterXSize();
    double adfGeoTransform[6];
    poSubDataset->GetGeoTransform(adfGeoTransform);
    if (adfGeoTransform[5] > 0) {
      adfGeoTransform[3] = adfGeoTransform[3] + adfGeoTransform[5] * rows;
      adfGeoTransform[5] = -adfGeoTransform[5];
    }
    // extent->PutCoord(adfGeoTransform[0], adfGeoTransform[3],
    //                  adfGeoTransform[0] + cols * adfGeoTransform[1],
    //                  adfGeoTransform[3] + adfGeoTransform[5] * rows);
    extent->MinX = adfGeoTransform[0];
    extent->MaxX = adfGeoTransform[0] + cols * adfGeoTransform[1];
    extent->MinY = adfGeoTransform[3] + rows * adfGeoTransform[5];
    extent->MaxY = adfGeoTransform[3];
    cellSize = fabs(adfGeoTransform[1]);
    poBand = poSubDataset->GetRasterBand(1);
  }
  *pVal = true;
  return;
}
void CGDALRasterReaderByPixel::GetPathName(std::string* pVal) { *pVal = this->lpszPathName; }
void CGDALRasterReaderByPixel::GetCurrentBand(PBand_T* pVal) { *pVal = pBand; }
void CGDALRasterReaderByPixel::GetRows(long* pVal) { *pVal = this->rows; }
void CGDALRasterReaderByPixel::GetCols(long* pVal) { *pVal = this->cols; }
void CGDALRasterReaderByPixel::GetExtent(OGREnvelope** pVal) { *pVal = this->extent; }
void CGDALRasterReaderByPixel::GetCellSize(double_t* pVal) { *pVal = this->cellSize; }
void CGDALRasterReaderByPixel::GetBandCount(long* pVal) {
  if (poDataset == nullptr) {
    *pVal = 0;
    return;
  } else if (metadata.size() == 0) {
    *pVal = poDataset->GetRasterCount();
  } else {
    if (poSubDataset == nullptr) {
      *pVal = 0;
      return;
    } else {
      *pVal = 1;
    }
  }
  return;
}
void CGDALRasterReaderByPixel::GetNoData(double_t* pVal) {
  if (poBand == nullptr) {
    *pVal = 0;
  } else {
    *pVal = poBand->GetNoDataValue();
  }
  return;
}
void CGDALRasterReaderByPixel::SetNoData(double_t pVal) {
  if (poBand != nullptr) {
    poBand->SetNoDataValue(pVal);
  }
}
void CGDALRasterReaderByPixel::GetMinMax(bool bApproxOK, double_t* min, double_t* max, bool* pVal) {
  if (poBand == nullptr) {
    *pVal = false;
    return;
  }

  double_t mean, stddev;
  CPLErr pErr = poBand->GetStatistics(bApproxOK, true, min, max, &mean, &stddev);
  if (pErr == CE_None)
    *pVal = true;
  else
    *pVal = false;
  return;
}
void CGDALRasterReaderByPixel::ComputeRasterMinMax(bool bApproxOK, double_t* min, double_t* max,
                                                   bool* pVal) {
  if (poBand == nullptr) {
    *pVal = false;
    return;
  }
  double dV[2];
  CPLErr pErr = poBand->ComputeRasterMinMax(bApproxOK, dV);
  if (pErr == CE_None) {
    *min = dV[0];
    *max = dV[1];
    *pVal = true;
  } else
    *pVal = false;
  return;
}
void CGDALRasterReaderByPixel::ComputeStatistics(bool bApproxOK, double_t* min, double_t* max,
                                                 double_t* mean, double_t* stddev, bool* pVal) {
  if (poBand == nullptr) {
    *pVal = false;
    return;
  }
  CPLErr pErr = poBand->ComputeStatistics(bApproxOK, min, max, mean, stddev, nullptr, nullptr);
  if (pErr == CE_None)
    *pVal = true;
  else
    *pVal = false;
  return;
}

void CGDALRasterReaderByPixel::GetStatistics(bool bApproxOK, double_t* min, double_t* max,
                                             double_t* mean, double_t* stddev, bool* pVal) {
  if (poBand == nullptr) {
    *pVal = false;
    return;
  }
  CPLErr pErr = poBand->GetStatistics(bApproxOK, true, min, max, mean, stddev);
  if (pErr == CE_None)
    *pVal = true;
  else
    *pVal = false;
  return;
}
// 获取参考系
void CGDALRasterReaderByPixel::GetSpatialReference(OGRSpatialReference** pVal) {
  OGRSpatialReference* pNew = nullptr;
  *pVal = pNew;
  GDALDataset* pDataset;
  if (poSubDataset != nullptr) {
    pDataset = poSubDataset;
  } else
    pDataset = poDataset;
  if (pDataset == nullptr) {
    return;
  }
  const char* info = pDataset->GetProjectionRef();
  if (info != nullptr) {
    pNew->importFromWkt(info);
  }
  return;
}
void CGDALRasterReaderByPixel::GetDataType(RasterDataType* pVal) {
  if (poBand == nullptr) {
    *pVal = RasterDataType::rdtUnknown;
    return;
  }
  *pVal = (RasterDataType)poBand->GetRasterDataType();
  return;
}
void CGDALRasterReaderByPixel::GetMetaData(std::vector<std::string>** pVal) {
  if (metadata.size() == 0) {
    *pVal = nullptr;
    return;
  }
  for (long k = 0; k < (int)metadata.size(); k++) {
    std::string sV = metadata.at(k);
    (*pVal)->push_back(sV);
  }
  return;
}

// TODO: 未实现
void CGDALRasterReaderByPixel::GetColorTable(std::vector<GDALColorTable>** pVal) {
  if (poBand == nullptr) {
    *pVal = nullptr;
    return;
  }
  GDALColorTable* gct = poBand->GetColorTable();
  if (gct == nullptr) {
    *pVal = nullptr;
    return;
  }
  int size = gct->GetColorEntryCount();
  // 遍历全部波段的颜色实体
  for (long k = 0; k < size; k++) {
    GDALColorEntry gce;
    gct->GetColorEntryAsRGB(k, &gce);
    // 拿到每个波段中的colorentry，保存到数组中
    GDALColorTable table;
    table.SetColorEntry(k, &gce);
    (*pVal)->push_back(table);
  }
  return;
}

/*-------------------CGDALRasterReaderByFileArray的实现--------------------------------*/
CGDALRasterReaderByFileArray::CGDALRasterReaderByFileArray() {
  rows = 0;
  cols = 0;
  CellSize = 0;
  poDataset = nullptr;
  poSubDataset = nullptr;
  poBand = nullptr;
  pData = nullptr;
  BufferSize = 0;
  poCT = nullptr;
  tpoCT = nullptr;
}

GDALColorTable* CGDALRasterReaderByFileArray::CreateColorTable() {
  GDALColorTable* pColorTable = nullptr;
  if (poBand == nullptr) {
    return nullptr;
  }
  GDALColorTable* gct = poBand->GetColorTable();
  if (gct == nullptr) return nullptr;
  pColorTable = new GDALColorTable(GPI_RGB);
  int Size = gct->GetColorEntryCount();
  for (int k = 0; k < Size; k++) {
    GDALColorEntry gce;
    gct->GetColorEntryAsRGB(k, &gce);
    pColorTable->SetColorEntry(k, &gce);
  }
  return pColorTable;
}

GDALRasterBand* CGDALRasterReaderByFileArray::GetGDALBand() { return this->poBand; }

// todo: 空实现
OGREnvelope CGDALRasterReaderByFileArray::TransformRect(OGRCoordinateTransformation* poCT,
                                                        OGREnvelope rt) {
  if (poCT == NULL) {
    return rt;
  }
  return rt;
}
OGREnvelope CGDALRasterReaderByFileArray::MapToPixelCoord(OGREnvelope MapExtent) {
  OGREnvelope PixelRect;
  // PixelRect.Left = (MapExtent.Left - XMin) / CellSize;
  // PixelRect.Right = (MapExtent.Right - XMin) / CellSize;
  // PixelRect.Top = (YMax - MapExtent.Top) / CellSize;
  // PixelRect.Bottom = (YMax - MapExtent.Bottom) / CellSize;
  return PixelRect;
}
OGREnvelope CGDALRasterReaderByFileArray::PixelToMapCoord(OGREnvelope PixelExtent) {
  OGREnvelope PaintExtent;
  // if (PixelExtent.Top > PixelExtent.Bottom) {
  //   double temp = PixelExtent.Top;
  //   PixelExtent.Top = PixelExtent.Bottom;
  //   PixelExtent.Bottom = temp;
  // }
  // PaintExtent.Left = PixelExtent.Left * CellSize + XMin;
  // PaintExtent.Right = (PixelExtent.Right + 1) * CellSize + XMin;
  // PaintExtent.Top = YMax - PixelExtent.Top * CellSize;
  // PaintExtent.Bottom = YMax - (PixelExtent.Bottom + 1) * CellSize;
  return PaintExtent;
}

void CGDALRasterReaderByFileArray::OpenRaster(std::string lpszPathName, bool* pVal) {
  // lpszPathName = (CString)PathName;
  if (poDataset != nullptr) {
    delete poDataset;
  }

  if (poSubDataset != nullptr) {
    delete poSubDataset;
  }
  poDataset = nullptr;
  poSubDataset = nullptr;
  rows = 0;
  cols = 0;
  XMin = XMax = YMin = YMax = 0;
  CellSize = 0;
  poBand = nullptr;
  BufferSize = 0;
  if (pData != nullptr) {
    // pData->Release();
    delete pData;
    pData = nullptr;
  }
  metadatas.clear();

  poDataset = (GDALDataset*)GDALOpen(lpszPathName.c_str(), GA_ReadOnly);
  if (poDataset == NULL) {
    *pVal = false;
    return;
  }
  const char* papszMetadata = GDALGetDriverShortName((GDALDriverH)poDataset);
  char** SUBDATASETS = GDALGetMetadata((GDALDatasetH)poDataset, "SUBDATASETS");
  if (CSLCount(SUBDATASETS) == 0) {
    rows = poDataset->GetRasterYSize();
    cols = poDataset->GetRasterXSize();
    double adfGeoTransform[6];
    poDataset->GetGeoTransform(adfGeoTransform);
    if (adfGeoTransform[5] > 0) {
      adfGeoTransform[3] = adfGeoTransform[3] + adfGeoTransform[5] * rows;
      adfGeoTransform[5] = -adfGeoTransform[5];
    }
    XMin = adfGeoTransform[0];
    YMin = adfGeoTransform[3] + adfGeoTransform[5] * rows;
    XMax = adfGeoTransform[0] + cols * adfGeoTransform[1];
    YMax = adfGeoTransform[3];
    CellSize = fabs(adfGeoTransform[1]);
  } else {
    for (int i = 0; SUBDATASETS[i] != nullptr; i++) {
      std::string tmpstr = SUBDATASETS[i];
      // 将每个SUBDATASETS按照等号分割
      int position = tmpstr.find("=", 0);
      // 取等号右边的属性值
      metadatas.push_back(tmpstr.substr(position));
    }
  }
  *pVal = true;
  return;
}

void CGDALRasterReaderByFileArray::PutRasterBand(PBand_T band, bool* pVal) {
  poBand = NULL;
  if (poDataset == NULL) {
    *pVal = false;
    return;
  }
  if (metadatas.size() == 0) {
    PBand_T vr = band;
    pBand = vr;
    poBand = poDataset->GetRasterBand(pBand.bandNumber);
  } else {
    rows = 0;
    cols = 0;
    CellSize = 0;
    XMin = XMax = YMin = YMax = 0;
    pBand = band;
    if (poSubDataset != NULL) delete poSubDataset;
    poSubDataset = NULL;
    poSubDataset = (GDALDataset*)GDALOpen(pBand.rasterPath, GA_ReadOnly);
    if (poSubDataset == NULL) {
      *pVal = false;
      return;
    }
    rows = poSubDataset->GetRasterYSize();
    cols = poSubDataset->GetRasterXSize();
    double adfGeoTransform[6];
    poSubDataset->GetGeoTransform(adfGeoTransform);
    if (adfGeoTransform[5] > 0) {
      adfGeoTransform[3] = adfGeoTransform[3] + adfGeoTransform[5] * rows;
      adfGeoTransform[5] = -adfGeoTransform[5];
    }
    XMin = adfGeoTransform[0];
    YMin = adfGeoTransform[3] + adfGeoTransform[5] * rows;
    XMax = adfGeoTransform[0] + cols * adfGeoTransform[1];
    YMax = adfGeoTransform[3];
    CellSize = fabs(adfGeoTransform[1]);
    poBand = poSubDataset->GetRasterBand(1);
  }
  *pVal = true;
  return;
}

void CGDALRasterReaderByFileArray::get_PathName(std::string* pVal) { *pVal = lpszPathName; }

void CGDALRasterReaderByFileArray::get_CurrentBand(PBand_T* pVal) { *pVal = pBand; }

void CGDALRasterReaderByFileArray::get_Rows(long* pVal) { *pVal = rows; }
void CGDALRasterReaderByFileArray::get_Cols(long* pVal) { *pVal = cols; }

// TODO: 返回了局部变量的引用，需要清理
void CGDALRasterReaderByFileArray::get_Extent(OGREnvelope** pVal) {
  OGREnvelope* pNew;
  pNew->MaxX = (*pVal)->MaxX;
  pNew->MaxY = (*pVal)->MaxY;
  pNew->MinX = (*pVal)->MinX;
  pNew->MinY = (*pVal)->MinY;
  *pVal = pNew;
  return;
}
void CGDALRasterReaderByFileArray::get_CellSize(double_t* pVal) { *pVal = this->CellSize; }
void CGDALRasterReaderByFileArray::get_BandCount(long* pVal) {
  if (poDataset == NULL) {
    *pVal = 0;
    return;
  } else if (metadatas.size() == 0) {
    *pVal = poDataset->GetRasterCount();
  } else {
    if (poSubDataset == NULL) {
      *pVal = 0;
      return;
    } else {
      *pVal = 1;
    }
  }
  return;
}
void CGDALRasterReaderByFileArray::get_NoData(double_t* pVal) {
  if (poBand == NULL) {
    *pVal = 0;
  } else {
    *pVal = poBand->GetNoDataValue();
  }
  return;
}

void CGDALRasterReaderByFileArray::put_NoData(double_t pVal) {
  if (poBand != NULL) poBand->SetNoDataValue(pVal);
  return;
}

void CGDALRasterReaderByFileArray::GetMinMax(bool bApproxOK, double_t* min, double_t* max,
                                             bool* pVal) {
  if (poBand == NULL) {
    *pVal = false;
    return;
  }
  double_t mean, stddev;
  CPLErr pErr = poBand->GetStatistics(bApproxOK, true, min, max, &mean, &stddev);
  if (pErr == CE_None)
    *pVal = true;
  else
    *pVal = false;
  return;
}

void CGDALRasterReaderByFileArray::GetStatistics(bool bApproxOK, double_t* min, double_t* max,
                                                 double_t* mean, double_t* stddev, bool* pVal) {
  if (poBand == NULL) {
    *pVal = false;
    return;
  }
  CPLErr pErr = poBand->GetStatistics(bApproxOK, true, min, max, mean, stddev);
  if (pErr == CE_None)
    *pVal = true;
  else
    *pVal = false;
  return;
}

void CGDALRasterReaderByFileArray::ComputeRasterMinMax(bool bApproxOK, double_t* min, double_t* max,
                                                       bool* pVal) {
  if (poBand == NULL) {
    *pVal = false;
    return;
  }
  double dV[2];
  CPLErr pErr = poBand->ComputeRasterMinMax(bApproxOK, dV);
  if (pErr == CE_None) {
    *min = dV[0];
    *max = dV[1];
    *pVal = true;
  } else
    *pVal = false;
  return;
}

void CGDALRasterReaderByFileArray::ComputeStatistics(bool bApproxOK, double_t* min, double_t* max,
                                                     double_t* mean, double_t* stddev, bool* pVal) {
  if (poBand == NULL) {
    *pVal = false;
    return;
  }
  CPLErr pErr = poBand->ComputeStatistics(bApproxOK, min, max, mean, stddev, NULL, NULL);
  if (pErr == CE_None)
    *pVal = true;
  else
    *pVal = false;
  return;
}

// TODO: 清理堆区资源
void CGDALRasterReaderByFileArray::GetSpatialReference(OGRSpatialReference** pVal) {
  GDALDataset* pDataset;
  if (poSubDataset != NULL) {
    pDataset = poSubDataset;
  } else {
    pDataset = poDataset;
  }
  if (pDataset == NULL) {
    return;
  }
  const char* info = pDataset->GetProjectionRef();
  if (info != NULL) {
    bool IsOk;
    // TODO: 有返回值
    (*pVal)->importFromWkt(info);
    return;
  }
}
void CGDALRasterReaderByFileArray::GetMetaData(std::vector<std::string>** pVal) {
  if (metadatas.size() == NULL) {
    *pVal = NULL;
    return;
  }
  std::vector<std::string> psa;
  long index;
  for (long k = 0; k < metadatas.size(); k++) {
    index = k;
    std::string sV = metadatas.at(k);
    psa.push_back(sV);
    // SafeArrayPutElement(psa, &index, sV.c_str());
  }
  *pVal = &psa;
  return;
}

void CGDALRasterReaderByFileArray::GetDataType(RasterDataType* pVal) {
  if (poBand == nullptr) {
    *pVal = RasterDataType::rdtUnknown;
    return;
  }
  *pVal = (RasterDataType)poBand->GetRasterDataType();
  return;
}

// TODO: 未实现GetColorTable
void CGDALRasterReaderByFileArray::GetColorTable(std::vector<GDALColorTable>** pVal) {
  if (poBand == nullptr) {
    *pVal = nullptr;
    return;
  }
  GDALColorTable* gct = poBand->GetColorTable();
  if (gct == nullptr) {
    *pVal = nullptr;
    return;
  }
  // int Size = gct->GetColorEntryCount();
  // // SAFEARRAY* psa;
  // std::vector<GDALColorTable> psa;
  // for (long k = 0; k < Size; k++) {
  //   GDALColorEntry gce;
  //   gct->GetColorEntryAsRGB(k, &gce);

  //   (*pVal)->push_back(gct);
  // }
  return;
}

void CGDALRasterReaderByFileArray::GetBlockData(long x1, long y1, long x2, long y2, long buffx,
                                                long buffy, IFileFloatArray** pVal) {
  if (poBand == NULL) {
    *pVal = NULL;
    return;
  }
  if (x1 > x2) {
    long temp = x1;
    x1 = x2;
    x2 = temp;
  }
  if (y1 > y2) {
    long temp = y1;
    y1 = y2;
    y2 = temp;
  }
  long buffer = buffx * buffy;
  if (buffer <= 0) {
    *pVal = nullptr;
    return;
  }
  if (buffer != BufferSize) {
    // if (pData != nullptr) {
    //   pData->Release();
    // }

    // TODO: 实例化一个pData类型的实例
    // this->pData = ;
    bool IsOk;
    pData->SetSize(buffer, poBand->GetNoDataValue(), &IsOk);
    if (!IsOk) {
      pData = nullptr;
      BufferSize = 0;
      *pVal = nullptr;
      return;
    }
    BufferSize = buffer;
  }
  float difx = (float)(x2 - x1 + 1) / buffx;
  float dify = (float)(y2 - y1 + 1) / buffy;
  float* fV = new float[x2 - x1 + 1];
  int CurrentRow = y1 - 1;
  long pos = 0;
  for (float row = y1 + dify / 2; row < y2 + 1; row += dify) {
    if (CurrentRow != (int)row) {
      poBand->RasterIO(GF_Read, x1, row, (x2 - x1 + 1), 1, fV, (x2 - x1 + 1), 1, GDT_Float32, 0, 0);
    }
    for (float col = x1 + difx / 2; col < x2 + 1; col += difx) {
      pData->SetValueAsFloat(pos, fV[(int)col]);
      pos++;
    }
    CurrentRow = row;
  }
  *pVal = pData;
  delete[] fV;
  return;
}

// TODO: 未实现GetInterpolatedBlockData
void CGDALRasterReaderByFileArray::GetInterpolatedBlockData(long x1, long y1, long x2, long y2,
                                                            long BuffX, long BuffY,
                                                            IFileFloatArray** pVal) {
  return;
}

void CGDALRasterReaderByFileArray::SetRefTargetSpatialReference(OGRSpatialReference* newVal) {
  if (poCT != nullptr) delete poCT;
  if (tpoCT != nullptr) delete tpoCT;
  poCT = nullptr;
  tpoCT = nullptr;
  if (newVal == nullptr) return;
  GDALDataset* pDataset;
  if (poSubDataset != nullptr)
    pDataset = poSubDataset;
  else
    pDataset = poDataset;
  if (pDataset == nullptr) {
    return;
  }
  const char* cinfo = pDataset->GetProjectionRef();
  OGRSpatialReference sp1;
  if (cinfo != nullptr) {
    sp1.importFromWkt(&cinfo);
  } else {
    return;
  }
  OGRSpatialReference sp2;
  char* wkt;
  newVal->exportToWkt(&wkt);
  if (wkt != nullptr) {
    sp2.importFromWkt(wkt);
  } else
    return;
  poCT = OGRCreateCoordinateTransformation(&sp2, &sp1);
  tpoCT = OGRCreateCoordinateTransformation(&sp1, &sp2);
  return;
}

void CGDALRasterReaderByFileArray::GetBlockDataByCoord(OGRPoint* LeftTop, double_t CellSize,
                                                       long Width, long Height, float NoData,
                                                       IFileFloatArray** pVal) {
  OGRPoint dpt;
  // LeftTop->GetCoord(&dpt.X, &dpt.Y);
  LeftTop->setX(dpt.getX());
  LeftTop->setY(dpt.getY());
  if (poCT == nullptr) {
    if (!GetDataBlock(dpt, CellSize, Width, Height, NoData))
      *pVal = nullptr;
    else
      *pVal = pData;
  } else {
    if (!GetDataBlockWithProj(dpt, CellSize, Width, Height, NoData))
      *pVal = nullptr;
    else
      *pVal = pData;
  }
  return;
}

void CGDALRasterReaderByFileArray::GetInterpolatedBlockDataByCoord(OGRPoint* LeftTop,
                                                                   double_t CellSize, long Width,
                                                                   long Height, float_t NoData,
                                                                   IFileFloatArray** pVal) {
  OGRPoint dpt;
  // LeftTop->GetCoord(&dpt.X, &dpt.Y);
  LeftTop->setX(dpt.getX());
  LeftTop->setY(dpt.getY());
  if (poCT == nullptr) {
    if (!GetInterpolatedDataBlock(dpt, CellSize, Width, Height, NoData))
      *pVal = nullptr;
    else
      *pVal = pData;
  } else {
    if (!GetInterpolatedDataBlockWithProj(dpt, CellSize, Width, Height, NoData))
      *pVal = nullptr;
    else
      *pVal = pData;
  }
  return;
}

// TODO: GetDataBlock未实现
bool CGDALRasterReaderByFileArray::GetDataBlock(long x1, long y1, long x2, long y2, long buffx,
                                                long buffy) {
  // DRect FullExtent = DRect(XMin, YMax, XMax, YMin);
  // float semiCellSize = cellSize / 2;
  // DRect CurrentExtent(LeftTop.X + semiCellSize, LeftTop.Y - semiCellSize,
  //                     LeftTop.X + cellSize * Width - semiCellSize,
  //                     LeftTop.Y - Height * cellSize + semiCellSize);
  // CRect ImageRect;
  // CRect TargetRect;
  // long nodata = poBand->GetNoDataValue();
  // LONG buffer = Width * Height;
  // if (buffer != BufferSize) {
  //   if (pData != NULL) {
  //     SafeArrayUnlock(pData);
  //     SafeArrayDestroy(pData);
  //   }
  //   BufferSize = buffer;
  //   pData = SafeArrayCreateVector(VT_R4, 0, BufferSize);
  //   SafeArrayLock(pData);
  // }
  // float* pvData = (float*)pData->pvData;
  // if (FullExtent.IsRectIn(CurrentExtent)) {
  //   ImageRect.left = (CurrentExtent.Left - XMin) / CellSize;
  //   ImageRect.top = (YMax - CurrentExtent.Top) / CellSize;
  //   ImageRect.right = (CurrentExtent.Right - XMin) / CellSize;
  //   ImageRect.bottom = (YMax - CurrentExtent.Bottom) / CellSize;
  //   if ((ImageRect.Width() + 1 >= Width) || (ImageRect.Height() + 1 >= Height)) {
  //     if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top, ImageRect.Width() + 1,
  //                          ImageRect.Height() + 1, pvData, Width, Height, GDT_Float32, 0, 0)
  //         != CE_None)
  //       return false;
  //     if (nodata != (long)NoData) {
  //       long Pos = 0;
  //       for (int i = 0; i < Height; i++) {
  //         for (int j = 0; j < Width; j++) {
  //           if ((long)pvData[Pos] == nodata) pvData[Pos] = NoData;
  //           Pos++;
  //         }
  //       }
  //     }
  //   } else {
  //     float* data = (float*)CPLMalloc(sizeof(float) * (ImageRect.Width() + 1) * sizeof(float)
  //                                     * (ImageRect.Height() + 1));
  //     if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top, ImageRect.Width() + 1,
  //                          ImageRect.Height() + 1, data, ImageRect.Width() + 1,
  //                          ImageRect.Height() + 1, GDT_Float32, 0, 0)
  //         != CE_None) {
  //       CPLFree(data);
  //       long Size = Width * Height;
  //       for (long k = 0; k < Size; k++) pvData[k] = NoData;
  //       return false;
  //     }
  //     DRect ext = PixelToMapCoord(ImageRect);
  //     float Y = LeftTop.Y - semiCellSize;
  //     long Pos = 0;
  //     float X;
  //     int posx, posy;
  //     long Posi;
  //     int W = ImageRect.Width();
  //     double rCellSize = CellSize;
  //     for (int i = 0; i < Height; i++) {
  //       X = LeftTop.X + semiCellSize;
  //       posy = (ext.Top - Y) / rCellSize;
  //       // if(posy>ImageRect.Height()) posy=ImageRect.Height();
  //       for (int j = 0; j < Width; j++) {
  //         posx = (X - ext.Left) / rCellSize;
  //         // if(posx>ImageRect.Width()) posx=ImageRect.Width();
  //         Posi = posy * (W + 1) + posx;
  //         if ((long)data[Posi] == nodata)
  //           pvData[Pos] = NoData;
  //         else
  //           pvData[Pos] = data[Posi];
  //         Pos++;
  //         X += cellSize;
  //       }
  //       Y -= cellSize;
  //     }
  //     CPLFree(data);
  //   }
  // } else {
  //   long Size = Width * Height;
  //   for (long k = 0; k < Size; k++) pvData[k] = NoData;
  //   if (!FullExtent.IntersectRect(CurrentExtent)) return true;
  //   CurrentExtent = FullExtent.Intersect(CurrentExtent);
  //   ImageRect.left = (CurrentExtent.Left - XMin) / CellSize;
  //   ImageRect.top = (YMax - CurrentExtent.Top) / CellSize;
  //   ImageRect.right = (CurrentExtent.Right - XMin) / CellSize;
  //   if (ImageRect.right >= cols) ImageRect.right = cols - 1;
  //   ImageRect.bottom = (YMax - CurrentExtent.Bottom) / CellSize;
  //   if (ImageRect.bottom >= rows) ImageRect.bottom = rows - 1;
  //   TargetRect.left = (CurrentExtent.Left - LeftTop.X) / cellSize;
  //   if (TargetRect.left < 0) TargetRect.left = 0;
  //   TargetRect.top = (LeftTop.Y - CurrentExtent.Top) / cellSize;
  //   if (TargetRect.top < 0) TargetRect.top = 0;
  //   TargetRect.right = (CurrentExtent.Right - LeftTop.X) / cellSize;
  //   if (TargetRect.right >= Width) TargetRect.right = Width - 1;
  //   TargetRect.bottom = (LeftTop.Y - CurrentExtent.Bottom) / cellSize;
  //   if (TargetRect.bottom >= Height) TargetRect.bottom = Height - 1;
  //   if ((ImageRect.Width() >= TargetRect.Width()) || (ImageRect.Height() >= TargetRect.Height()))
  //   {
  //     float* data = (float*)CPLMalloc(sizeof(float) * (TargetRect.Width() + 1) * sizeof(float)
  //                                     * (TargetRect.Height() + 1));
  //     if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top, ImageRect.Width() + 1,
  //                          ImageRect.Height() + 1, data, TargetRect.Width() + 1,
  //                          TargetRect.Height() + 1, GDT_Float32, 0, 0)
  //         != CE_None) {
  //       CPLFree(data);
  //       return false;
  //     }
  //     long Pos;
  //     long Posi = 0;
  //     for (int i = TargetRect.top; i <= TargetRect.bottom; i++) {
  //       Pos = i * Width + TargetRect.left;
  //       for (int j = TargetRect.left; j <= TargetRect.right; j++) {
  //         if ((long)data[Posi] == nodata)
  //           pvData[Pos] = NoData;
  //         else
  //           pvData[Pos] = data[Posi];
  //         Pos++;
  //         Posi++;
  //       }
  //     }
  //     CPLFree(data);
  //   } else {
  //     float* data = (float*)CPLMalloc(sizeof(float) * (ImageRect.Width() + 1) * sizeof(float)
  //                                     * (ImageRect.Height() + 1));
  //     if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top, ImageRect.Width() + 1,
  //                          ImageRect.Height() + 1, data, ImageRect.Width() + 1,
  //                          ImageRect.Height() + 1, GDT_Float32, 0, 0)
  //         != CE_None) {
  //       CPLFree(data);
  //       return false;
  //     }
  //     long Pos;
  //     DRect ext = PixelToMapCoord(ImageRect);
  //     float Y = LeftTop.Y - semiCellSize - TargetRect.top * cellSize;
  //     float X;
  //     int posx, posy;
  //     long Posi;
  //     int W = ImageRect.Width();
  //     double rCellSize = CellSize;
  //     for (int i = TargetRect.top; i <= TargetRect.bottom; i++) {
  //       X = LeftTop.X + semiCellSize + TargetRect.left * cellSize;
  //       Pos = i * Width + TargetRect.left;
  //       posy = (ext.Top - Y) / rCellSize;
  //       if (posy > ImageRect.Height()) posy = ImageRect.Height();
  //       for (int j = TargetRect.left; j <= TargetRect.right; j++) {
  //         posx = (X - ext.Left) / rCellSize;
  //         if (posx > W) posx = W;
  //         Posi = posy * (W + 1) + posx;
  //         if ((long)data[Posi] == nodata)
  //           pvData[Pos] = NoData;
  //         else
  //           pvData[Pos] = data[Posi];
  //         Pos++;
  //         X += cellSize;
  //       }
  //       Y -= cellSize;
  //     }
  //     CPLFree(data);
  //   }
  // }
  // return true;
  return true;
}

bool CGDALRasterReaderByFileArray::ReadDataBlock(long x1, long y1, long x2, long y2, long buffx,
                                                 long buffy, IFileFloatArray* pData) {
  if (poBand == nullptr) {
    return false;
  }
  if (x1 > x2) {
    long temp = x1;
    x1 = x2;
    x2 = temp;
  }
  if (y1 > y2) {
    long temp = y1;
    y1 = y2;
    y2 = temp;
  }
  long buffer = buffx * buffy;
  if (buffer <= 0) {
    return false;
  }
  float difx = (float)(x2 - x1 + 1) / buffx;
  float dify = (float)(y2 - y1 + 1) / buffy;
  float* fV = new float[x2 - x1 + 1];
  int CurrentRow = y1 - 1;
  long pos = 0;
  for (float row = y1 + dify / 2; row < y2 + 1; row += dify) {
    if (CurrentRow != (int)row) {
      poBand->RasterIO(GF_Read, x1, row, (x2 - x1 + 1), 1, fV, (x2 - x1 + 1), 1, GDT_Float32, 0, 0);
    }
    for (float col = x1 + difx / 2; col < x2 + 1; col += difx) {
      pData->SetValueAsFloat(pos, fV[(int)col - x1]);
      pos++;
    }
    CurrentRow = row;
  }
  delete[] fV;
  return true;
}

IFileFloatArray* CGDALRasterReaderByFileArray::GetDataBlock(long x1, long y1, long x2, long y2,
                                                            long buffx, long buffy) {
  if (poBand == NULL) {
    return NULL;
  }
  if (x1 > x2) {
    long temp = x1;
    x1 = x2;
    x2 = temp;
  }
  if (y1 > y2) {
    long temp = y1;
    y1 = y2;
    y2 = temp;
  }
  long buffer = buffx * buffy;
  if (buffer <= 0) {
    return NULL;
  }
  // pData要返回，需要在堆上新建
  IFileFloatArray* pData;
  bool IsOk;
  pData->SetSize(buffer, poBand->GetNoDataValue(), &IsOk);
  if (!IsOk) {
    pData = NULL;
    return NULL;
  }
  float difx = (float)(x2 - x1 + 1) / buffx;
  float dify = (float)(y2 - y1 + 1) / buffy;
  float* fV = new float[x2 - x1 + 1];
  int CurrentRow = y1 - 1;
  long pos = 0;
  for (float row = y1 + dify / 2; row < y2 + 1; row += dify) {
    if (CurrentRow != (int)row) {
      poBand->RasterIO(GF_Read, x1, row, (x2 - x1 + 1), 1, fV, (x2 - x1 + 1), 1, GDT_Float32, 0, 0);
    }
    for (float col = x1 + difx / 2; col < x2 + 1; col += difx) {
      pData->SetValueAsFloat(pos, fV[(int)col - x1]);
      pos++;
    }
    CurrentRow = row;
  }
  delete[] fV;
  return pData;
}

bool CGDALRasterReaderByFileArray::GetDataBlock(OGRPoint LeftTop, float cellSize, int Width,
                                                int Height, float NoData) {
  // DRect FullExtent = DRect(XMin, YMax, XMax, YMin);
  OGREnvelope FullExtent;
  FullExtent.MinX = XMin;
  FullExtent.MaxY = YMax;
  FullExtent.MaxX = XMax;
  FullExtent.MinY = YMin;
  float semiCellSize = cellSize / 2;

  // DRect CurrentExtent(LeftTop.getX() + semiCellSize, LeftTop.getY() - semiCellSize,
  //                     LeftTop.getX() + cellSize * Width - semiCellSize,
  //                     LeftTop.getY() - Height * cellSize + semiCellSize);

  // 左/上/右/下的顺序
  OGREnvelope CurrentExtent;
  CurrentExtent.MinX = LeftTop.getX() + semiCellSize;
  CurrentExtent.MaxY = LeftTop.getY() - semiCellSize;
  CurrentExtent.MaxX = LeftTop.getX() + cellSize * Width - semiCellSize;
  CurrentExtent.MinY = LeftTop.getY() - Height * cellSize + semiCellSize;

  CRect ImageRect;
  CRect TargetRect;
  long nodata = poBand->GetNoDataValue();
  long buffer = Width * Height;
  if (buffer != BufferSize) {
    if (pData != NULL) {
      pData = NULL;
    }

    // TODO: 实例化pData
    bool IsOk;
    pData->SetSize(buffer, NoData, &IsOk);
    if (!IsOk) {
      pData = NULL;
      BufferSize = 0;
      return;
    }
    BufferSize = buffer;
  }

  if (FullExtent.Contains(CurrentExtent)) {
    ImageRect.left = (CurrentExtent.MinX - XMin) / CellSize;
    ImageRect.top = (YMax - CurrentExtent.Top) / CellSize;
    ImageRect.right = (CurrentExtent.Right - XMin) / CellSize;
    if (ImageRect.right >= cols) ImageRect.right = cols - 1;
    ImageRect.bottom = (YMax - CurrentExtent.Bottom) / CellSize;
    if (ImageRect.bottom >= rows) ImageRect.bottom = rows - 1;
    if ((ImageRect.Width() + 1 >= Width) || (ImageRect.Height() + 1 >= Height)) {
      if (!ReadDataBlock(ImageRect.left, ImageRect.top, ImageRect.right, ImageRect.bottom, Width,
                         Height, pData))
        return false;
      if (nodata != (long)NoData) {
        long Pos = 0;
        for (int i = 0; i < Height; i++) {
          for (int j = 0; j < Width; j++) {
            float v;
            pData->GetValueAsFloat(Pos, &v);
            if ((long)v == nodata) pData->SetValueAsFloat(Pos, NoData);
            Pos++;
          }
          // if (progress != NULL) progress->SetPos((float)i / Height * 100);
        }
      }
    } else {
      IFileFloatArray* data
          = GetDataBlock(ImageRect.left, ImageRect.top, ImageRect.right, ImageRect.bottom,
                         ImageRect.Width() + 1, ImageRect.Height() + 1);
      if (data == NULL) {
        pData->SetDefaultValue(NoData);
        return false;
      }
      DRect ext = PixelToMapCoord(ImageRect);
      float Y = LeftTop.getY() - semiCellSize;
      long Pos = 0;
      float X;
      int posx, posy;
      long Posi;
      int W = ImageRect.Width();
      double rCellSize = CellSize;
      for (int i = 0; i < Height; i++) {
        X = LeftTop.X + semiCellSize;
        posy = (ext.Top - Y) / rCellSize;
        // if(posy>ImageRect.Height()) posy=ImageRect.Height();
        for (int j = 0; j < Width; j++) {
          posx = (X - ext.Left) / rCellSize;
          // if(posx>ImageRect.Width()) posx=ImageRect.Width();
          Posi = posy * (W + 1) + posx;
          float v;
          data->GetValueAsFloat(Posi, &v);
          if ((long)v == nodata)
            pData->SetValueAsFloat(Pos, NoData);
          else
            pData->SetValueAsFloat(Pos, v);
          Pos++;
          X += cellSize;
        }

        Y -= cellSize;
      }
    }
  } else {
    pData->SetDefaultValue(NoData);
    long Size = Width * Height;
    if (!FullExtent.IntersectRect(CurrentExtent)) return true;
    CurrentExtent = FullExtent.Intersect(CurrentExtent);
    ImageRect.left = (CurrentExtent.Left - XMin) / CellSize;
    ImageRect.top = (YMax - CurrentExtent.Top) / CellSize;
    ImageRect.right = (CurrentExtent.Right - XMin) / CellSize;
    if (ImageRect.right >= cols) ImageRect.right = cols - 1;
    ImageRect.bottom = (YMax - CurrentExtent.Bottom) / CellSize;
    if (ImageRect.bottom >= rows) ImageRect.bottom = rows - 1;
    TargetRect.left = (CurrentExtent.Left - LeftTop.X) / cellSize;
    if (TargetRect.left < 0) TargetRect.left = 0;
    TargetRect.top = (LeftTop.Y - CurrentExtent.Top) / cellSize;
    if (TargetRect.top < 0) TargetRect.top = 0;
    TargetRect.right = (CurrentExtent.Right - LeftTop.X) / cellSize;
    if (TargetRect.right >= Width) TargetRect.right = Width - 1;
    TargetRect.bottom = (LeftTop.Y - CurrentExtent.Bottom) / cellSize;
    if (TargetRect.bottom >= Height) TargetRect.bottom = Height - 1;
    if ((ImageRect.Width() >= TargetRect.Width()) || (ImageRect.Height() >= TargetRect.Height())) {
      IFileFloatArray* data
          = GetDataBlock(ImageRect.left, ImageRect.top, ImageRect.right, ImageRect.bottom,
                         TargetRect.Width() + 1, TargetRect.Height() + 1, progress);
      if (data == NULL) return false;
      long Pos;
      long Posi = 0;
      for (int i = TargetRect.top; i <= TargetRect.bottom; i++) {
        Pos = i * Width + TargetRect.left;
        for (int j = TargetRect.left; j <= TargetRect.right; j++) {
          float v;
          data->GetValueAsFloat(Posi, &v);
          if ((long)v == nodata)
            pData->SetValueAsFloat(Pos, NoData);
          else
            pData->SetValueAsFloat(Pos, v);
          Pos++;
          Posi++;
        }
        if (progress != NULL)
          progress->SetPos((float)(i - TargetRect.top) / (TargetRect.bottom - TargetRect.top)
                           * 100);
      }
      // data->Release();
    } else {
      IFileFloatArray* data
          = GetDataBlock(ImageRect.left, ImageRect.top, ImageRect.right, ImageRect.bottom,
                         ImageRect.Width() + 1, ImageRect.Height() + 1);
      if (data == NULL) return false;
      long Pos;
      DRect ext = PixelToMapCoord(ImageRect);
      float Y = LeftTop.Y - semiCellSize - TargetRect.top * cellSize;
      float X;
      int posx, posy;
      long Posi;
      int W = ImageRect.Width();
      double rCellSize = CellSize;
      for (int i = TargetRect.top; i <= TargetRect.bottom; i++) {
        X = LeftTop.X + semiCellSize + TargetRect.left * cellSize;
        Pos = i * Width + TargetRect.left;
        posy = (ext.Top - Y) / rCellSize;
        if (posy > ImageRect.Height()) posy = ImageRect.Height();
        for (int j = TargetRect.left; j <= TargetRect.right; j++) {
          posx = (X - ext.Left) / rCellSize;
          if (posx > W) posx = W;
          Posi = posy * (W + 1) + posx;
          float v;
          data->GetValueAsFloat(Posi, &v);
          if ((long)v == nodata)
            pData->SetValueAsFloat(Pos, NoData);
          else
            pData->SetValueAsFloat(Pos, v);
          Pos++;
          X += cellSize;
        }
        Y -= cellSize;
        if (progress != NULL)
          progress->SetPos((float)(i - TargetRect.top) / (TargetRect.bottom - TargetRect.top)
                           * 100);
      }
      // data->Release();
    }
  }
  return true;
}