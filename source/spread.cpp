#include <fmt/format.h>
#include <spread/spread.h>

#include <iostream>
#include <vector>

using namespace spread;

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
  // ::CoCreateInstance(CLSID_PointDef, NULL, CLSCTX_INPROC_SERVER, IID_IPoint, (void**)&ppt);
  double_t left, top;
  // Extent->GetLeft(&left);
  // Extent->GetTop(&top);
  left = Extent->MinY;
  top = Extent->MaxX;
  // ppt->PutCoord(left, top);
  ppt->setX(top);
  ppt->setY(left);
  pPro->GetCellSize(&cellSize);
  if (pEnvi == NULL) {
    // ::CoCreateInstance(CLSID_AnalyseEnvi, NULL, CLSCTX_INPROC_SERVER, IID_IAnalyseEnvi,
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
    if (lt != NULL) {
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
  pElevData = NULL;
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
    if (otherDatas.at(k) != NULL) {
      std::vector<IFileFloatArray*>::iterator pos = otherDatas.begin() + k;
      otherDatas.erase(pos);
    }
    otherDatas.push_back(pArray);
    if (pArray == NULL) {
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
  GDALColorTable* pColorTable = NULL;
  if (poBand == NULL) {
    return NULL;
  }
  GDALColorTable* gct = poBand->GetColorTable();
  if (gct == NULL) {
    return NULL;
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
  if (poDataset != NULL) {
    delete poDataset;
  }
  if (poSubDataset != NULL) {
    delete poSubDataset;
  }
  poDataset = NULL;
  poSubDataset = NULL;
  rows = 0;
  cols = 0;
  cellSize = 0;
  poBand = NULL;
  // extent->PutCoord(0, 0, 0, 0);
  extent->MaxX = 0;
  extent->MaxY = 0;
  extent->MinX = 0;
  extent->MinY = 0;
  metadata.clear();

  // 打开栅格数据文件，获得Dataset对象
  poDataset = (GDALDataset*)GDALOpen(lpszPathName.c_str(), GA_ReadOnly);
  if (poDataset == NULL) {
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
    for (int i = 0; SUBDATASETS[i] != NULL; i++) {
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
  if (poBand == NULL) {
    *data = -INT_MAX;
    return;
  }
  CPLErr errorNumber = poBand->RasterIO(GF_Read, col, row, 1, 1, data, 1, 1, GDT_Float32, 0, 0);
  std::cout << "获取结果标志：" << errorNumber << std::endl;
}
void CGDALRasterReaderByPixel::SetRasterBand(PBand_T pBand, bool* pVal) {
  poBand = NULL;
  if (poDataset == NULL) {
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

    if (poSubDataset != NULL) {
      delete poSubDataset;
    }
    poSubDataset = NULL;
    // _bstr_t path = vr.bstrVal;
    // char* cpath = path;

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
  if (poDataset == NULL) {
    *pVal = 0;
    return;
  } else if (metadata.size() == 0) {
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
void CGDALRasterReaderByPixel::GetNoData(double_t* pVal) {
  if (poBand == NULL) {
    *pVal = 0;
  } else {
    *pVal = poBand->GetNoDataValue();
  }
  return;
}
void CGDALRasterReaderByPixel::SetNoData(double_t pVal) {
  if (poBand != NULL) {
    poBand->SetNoDataValue(pVal);
  }
}
void CGDALRasterReaderByPixel::GetMinMax(bool bApproxOK, double_t* min, double_t* max, bool* pVal) {
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
void CGDALRasterReaderByPixel::ComputeRasterMinMax(bool bApproxOK, double_t* min, double_t* max,
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
void CGDALRasterReaderByPixel::ComputeStatistics(bool bApproxOK, double_t* min, double_t* max,
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
// 获取参考系
void CGDALRasterReaderByPixel::GetSpatialReference(OGRSpatialReference** pVal) {
  OGRSpatialReference* pNew = nullptr;
  *pVal = pNew;
  GDALDataset* pDataset;
  if (poSubDataset != NULL) {
    pDataset = poSubDataset;
  } else
    pDataset = poDataset;
  if (pDataset == NULL) {
    return;
  }
  const char* info = pDataset->GetProjectionRef();
  if (info != NULL) {
    pNew->importFromWkt(info);
  }
  return;
}
void CGDALRasterReaderByPixel::GetDataType(RasterDataType* pVal) {
  if (poBand == NULL) {
    *pVal = RasterDataType::rdtUnknown;
    return;
  }
  *pVal = (RasterDataType)poBand->GetRasterDataType();
  return;
}
void CGDALRasterReaderByPixel::GetMetaData(std::vector<std::string>** pVal) {
  if (metadata.size() == 0) {
    *pVal = NULL;
    return;
  }
  for (long k = 0; k < (int)metadata.size(); k++) {
    std::string sV = metadata.at(k);
    (*pVal)->push_back(sV);
  }
  return;
}
void CGDALRasterReaderByPixel::GetColorTable(std::vector<GDALColorTable>** pVal) {
  if (poBand == NULL) {
    *pVal = NULL;
    return;
  }
  GDALColorTable* gct = poBand->GetColorTable();
  if (gct == NULL) {
    *pVal = NULL;
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