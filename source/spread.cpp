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
  if (!pElevs) {
    this->pElevs = new CGDALRasterReaderByFileArray;
  }

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

  IGDALRasterProperties* pPro = dynamic_cast<CGDALRasterReaderByFileArray*>(this->pElevs);
  // pElevs->QueryInterface(IID_IGDALRasterProperties, (void**)&pPro);
  pPro->GetNoData(&noData);

  pPro->GetCols(&cols);
  pPro->GetRows(&rows);
  OGREnvelope* Extent;
  pPro->GetExtent(&Extent);
  OGRPoint* ppt = new OGRPoint;
  // ::CoCreateInstance(CLSID_PointDef, nullptr, CLSCTX_INPROC_SERVER, IID_IPoint, (void**)&ppt);
  double_t left, top;
  // Extent->GetLeft(&left);
  // Extent->GetTop(&top);
  left = Extent->MinX;
  top = Extent->MaxY;
  // ppt->PutCoord(left, top);
  ppt->setX(left);
  ppt->setY(top);
  pPro->GetCellSize(&cellSize);
  if (pEnvi == nullptr) {
    // ::CoCreateInstance(CLSID_AnalyseEnvi, nullptr, CLSCTX_INPROC_SERVER, IID_IAnalyseEnvi,
    //                    (LPVOID*)&pEnvi);  // 得到接口指针
    pEnvi = new AnalyseEnvironment;
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
  // TODO: 在堆中分配空间
  pElevData = new IFileFloatArray;
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
int CRect::Width() { return this->right - this->left; }
int CRect::Height() { return this->bottom - this->top; }

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
  pStations = new CStations;
}
bool CFieldStrengthAnalyse::FieldStrengthAnalyse(std::string savePath, RasterCreateFileType type) {
  errorInfo = "";
  long Count = 0;
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
  pData = new IFileFloatArray;
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
  // delete pData;
  return true;
}

void CFieldStrengthAnalyse::ComputeOneStation(Station& para) {
  std::cout << "计算" + GetModelName() + "场强" << std::endl;

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
  leftTop = new OGRPoint;
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
void CStations::AddStation(Station* station) { this->stations.push_back(station); }
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
  PixelRect.MinX = (MapExtent.MinX - XMin) / CellSize;
  PixelRect.MinY = (MapExtent.MaxX - XMin) / CellSize;
  PixelRect.MaxY = (YMax - MapExtent.MaxY) / CellSize;
  PixelRect.MinY = (YMax - MapExtent.MinY) / CellSize;
  return PixelRect;
}
OGREnvelope CGDALRasterReaderByFileArray::PixelToMapCoord(OGREnvelope PixelExtent) {
  OGREnvelope PaintExtent;
  if (PixelExtent.MaxY > PixelExtent.MinY) {
    double_t temp = PixelExtent.MaxY;
    PixelExtent.MaxY = PixelExtent.MinY;
    PixelExtent.MinY = temp;
  }
  PaintExtent.MinX = PixelExtent.MinX * CellSize + XMin;
  PaintExtent.MaxX = (PixelExtent.MaxX + 1) * CellSize + XMin;
  PaintExtent.MaxY = YMax - PixelExtent.MaxY * CellSize;
  PaintExtent.MinY = YMax - (PixelExtent.MinY + 1) * CellSize;
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
  GDALAllRegister();
  poDataset = (GDALDataset*)GDALOpen(lpszPathName.c_str(), GA_ReadOnly);
  if (poDataset == NULL) {
    *pVal = false;
    return;
  }
  // const char* papszMetadata = GDALGetDriverShortName((GDALDriverH)poDataset);
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

void CGDALRasterReaderByFileArray::SetRasterBand(PBand_T band, bool* pVal) {
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

void CGDALRasterReaderByFileArray::GetPathName(std::string* pVal) { *pVal = lpszPathName; }

void CGDALRasterReaderByFileArray::GetCurrentBand(PBand_T* pVal) { *pVal = pBand; }

void CGDALRasterReaderByFileArray::GetRows(long* pVal) { *pVal = rows; }
void CGDALRasterReaderByFileArray::GetCols(long* pVal) { *pVal = cols; }

// TODO: 返回了局部变量的引用，需要清理
void CGDALRasterReaderByFileArray::GetExtent(OGREnvelope** pVal) {
  OGREnvelope* pNew = new OGREnvelope();
  pNew->MaxX = XMax;
  pNew->MaxY = YMax;
  pNew->MinX = XMin;
  pNew->MinY = YMin;
  *pVal = pNew;
  return;
}
void CGDALRasterReaderByFileArray::GetCellSize(double_t* pVal) { *pVal = this->CellSize; }
void CGDALRasterReaderByFileArray::GetBandCount(long* pVal) {
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
void CGDALRasterReaderByFileArray::GetNoData(double_t* pVal) {
  if (poBand == NULL) {
    *pVal = 0;
  } else {
    *pVal = poBand->GetNoDataValue();
  }
  return;
}

void CGDALRasterReaderByFileArray::SetNoData(double_t pVal) {
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
    int IsOk;
    // TODO: 有返回值
    IsOk = (*pVal)->importFromWkt(info);
    if (IsOk != 0) {
      return;
    }
    return;
  }
}
void CGDALRasterReaderByFileArray::GetMetaData(std::vector<std::string>** pVal) {
  if (metadatas.size() == 0L) {
    *pVal = NULL;
    return;
  }
  std::vector<std::string> psa;
  // long index;
  for (long k = 0; k < (long)metadatas.size(); k++) {
    // index = k;
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
      CPLErr err = poBand->RasterIO(GF_Read, x1, row, (x2 - x1 + 1), 1, fV, (x2 - x1 + 1), 1,
                                    GDT_Float32, 0, 0);
      if (err != CE_None) {
        return;
      }
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
  // AFX_MANAGE_STATE(AfxGetStaticModuleState());
  int Width = BuffX;
  int Height = BuffY;
  CRect ImageRect;
  long nodata = poBand->GetNoDataValue();
  long buffer = Width * Height;
  if (buffer != BufferSize) {
    if (pData != NULL) {
      // pData->Release();
      pData = NULL;
    }
    // ::CoCreateInstance(CLSID_FileFloatArray, NULL, CLSCTX_INPROC_SERVER, IID_IFileFloatArray,
    //                    (LPVOID*)&pData);
    bool IsOk;
    pData->SetSize(buffer, nodata, &IsOk);
    if (!IsOk) {
      // pData->Release();
      pData = NULL;
      BufferSize = 0;
      return;
    }
    BufferSize = buffer;
  }
  ImageRect = CRect(x1, y1, x2, y2);
  // ImageRect.NormalizeRect();
  int nTemp;
  if (ImageRect.left > ImageRect.right) {
    nTemp = ImageRect.left;
    ImageRect.left = ImageRect.right;
    ImageRect.right = nTemp;
  }
  if (ImageRect.top > ImageRect.bottom) {
    nTemp = ImageRect.top;
    ImageRect.top = ImageRect.bottom;
    ImageRect.bottom = nTemp;
  }

  x1 = ImageRect.left;
  x2 = ImageRect.right;
  y1 = ImageRect.top;
  y2 = ImageRect.bottom;
  if ((ImageRect.Width() + 1 >= Width) || (ImageRect.Height() + 1 >= Height)) {
    if (!ReadDataBlock(ImageRect.left, ImageRect.top, ImageRect.right, ImageRect.bottom, Width,
                       Height, pData)) {
      *pVal = NULL;
      return;
    }
  } else {
    // ImageRect.InflateRect(1, 1, 1, 1);
    ImageRect.left -= 1;
    ImageRect.top -= 1;
    ImageRect.right += 1;
    ImageRect.bottom += 1;
    if (ImageRect.left < 0) ImageRect.left++;
    if (ImageRect.right >= cols) ImageRect.right--;
    if (ImageRect.top < 0) ImageRect.top++;
    if (ImageRect.bottom >= rows) ImageRect.bottom--;
    IFileFloatArray* data
        = GetDataBlock(ImageRect.left, ImageRect.top, ImageRect.right, ImageRect.bottom,
                       ImageRect.Width() + 1, ImageRect.Height() + 1);
    if (data == NULL) {
      pData->SetDefaultValue(nodata);
      return;
    }
    long Pos = 0;
    float posx, posy;
    int iposx, iposy, iposy1;
    int State;  // 0--None;1--Left;2--Right;
    long Posi;
    int W = ImageRect.Width();
    int H = ImageRect.Height();
    float ix, ix2;
    float ratiox = (float)(W + 1) / Width;
    float ratioy = (float)(H + 1) / Height;
    // if (progress != NULL) progress->BeginProgress(CComBSTR("拷贝数据"));
    for (int i = 0; i < Height; i++) {
      posy = y1 + ratioy * (0.5 + i);
      iposy = posy;
      if (posy - iposy < 0.5) {
        iposy1 = iposy - 1;
        if (iposy1 < 0) iposy1 = iposy;
      } else {
        iposy1 = iposy + 1;
        if (iposy1 > H) iposy1 = iposy;
      }
      for (int j = 0; j < Width; j++) {
        posx = x1 + ratiox * (0.5 + j);
        iposx = posx;
        Posi = iposy * (W + 1) + iposx;
        float v0;
        data->GetValueAsFloat(Posi, &v0);
        if ((long)v0 == nodata)
          ix = nodata;
        else {
          if (posx - iposx < 0.5) {
            if (iposx > 0) {
              float v;
              data->GetValueAsFloat(Posi - 1, &v);
              if ((long)v != nodata)
                State = 1;
              else
                State = 2;
            } else
              State = 2;
            if (State == 2) {
              if (iposx + 1 > W)
                State = 0;
              else {
                float v;
                data->GetValueAsFloat(Posi + 1, &v);
                if ((long)v == nodata) State = 0;
              }
            }
          } else {
            if (iposx + 1 <= W) {
              float v;
              data->GetValueAsFloat(Posi + 1, &v);
              if ((long)v != nodata)
                State = 2;
              else
                State = 1;
            } else
              State = 1;
            if (State == 1) {
              if (iposx < 1)
                State = 0;
              else {
                float v;
                data->GetValueAsFloat(Posi - 1, &v);
                if ((long)v == nodata) State = 0;
              }
            }
          }
          switch (State) {
            case 0: {
              ix = v0;
              break;
            }
            case 1: {
              float v;
              data->GetValueAsFloat(Posi - 1, &v);
              ix = (posx - iposx + 0.5) * (v0 - v) + v;
              break;
            }
            case 2: {
              float v;
              data->GetValueAsFloat(Posi + 1, &v);
              ix = (posx - iposx - 0.5) * (v - v0) + v0;
              break;
            }
          }
        }
        if (iposy1 == iposy) {
          pData->SetValueAsFloat(Pos, ix);
          Pos++;
          continue;
        }
        Posi = iposy1 * (W + 1) + iposx;
        data->GetValueAsFloat(Posi, &v0);
        if ((long)v0 == nodata)
          ix2 = nodata;
        else {
          if (posx - iposx < 0.5) {
            if (iposx > 0) {
              float v;
              data->GetValueAsFloat(Posi - 1, &v);
              if ((long)v != nodata)
                State = 1;
              else
                State = 2;
            } else
              State = 2;
            if (State == 2) {
              if (iposx + 1 > W)
                State = 0;
              else {
                float v;
                data->GetValueAsFloat(Posi + 1, &v);
                if ((long)v == nodata) State = 0;
              }
            }
          } else {
            if (iposx + 1 <= W) {
              float v;
              data->GetValueAsFloat(Posi + 1, &v);
              if ((long)v != nodata)
                State = 2;
              else
                State = 1;
            } else
              State = 1;
            if (State == 1) {
              if (iposx < 1)
                State = 0;
              else {
                float v;
                data->GetValueAsFloat(Posi - 1, &v);
                if ((long)v == nodata) State = 0;
              }
            }
          }
          switch (State) {
            case 0: {
              ix2 = v0;
              break;
            }
            case 1: {
              float v;
              data->GetValueAsFloat(Posi - 1, &v);
              ix2 = (posx - iposx + 0.5) * (v0 - v) + v;
              break;
            }
            case 2: {
              float v;
              data->GetValueAsFloat(Posi + 1, &v);
              ix2 = (posx - iposx - 0.5) * (v - v0) + v0;
              break;
            }
          }
        }
        State = 0;
        if ((long)ix == nodata) State = 1;
        if ((long)ix2 == nodata) State += 2;
        switch (State) {
          case 0:
            pData->SetValueAsFloat(Pos, (ix2 - ix) * (posy - iposy - 0.5) / (iposy1 - iposy) + ix);
            break;
          case 1:
            pData->SetValueAsFloat(Pos, ix2);
            break;
          case 2:
            pData->SetValueAsFloat(Pos, ix);
            break;
          case 3:
            pData->SetValueAsFloat(Pos, nodata);
            break;
        }
        Pos++;
      }
      // if (progress != NULL) progress->SetPos((float)i / Height * 100);
    }
    // data->Release();
  }
  *pVal = pData;
  // if (*pVal != NULL) (*pVal)->AddRef();
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
  // LeftTop->setX(dpt.getX());
  dpt.setX(LeftTop->getX());
  dpt.setY(LeftTop->getY());
  // LeftTop->setY(dpt.getY());
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

bool CGDALRasterReaderByFileArray::GetDataBlock(OGRPoint LeftTop, float cellSize, int Width,
                                                int Height, float NoData) {
  OGREnvelope FullExtent;
  FullExtent.MinX = XMin;
  FullExtent.MinY = YMin;
  FullExtent.MaxX = XMax;
  FullExtent.MaxY = YMax;
  // DRect(XMin, YMax, XMax, YMin);
  float semiCellSize = cellSize / 2;
  OGREnvelope CurrentExtent;
  CurrentExtent.MinX = LeftTop.getX() + semiCellSize;
  CurrentExtent.MinY = LeftTop.getY() - Height * cellSize + semiCellSize;
  CurrentExtent.MaxX = LeftTop.getX() + cellSize * Width - semiCellSize;
  CurrentExtent.MaxY = LeftTop.getY() - semiCellSize;
  // (LeftTop.X + semiCellSize, LeftTop.Y - semiCellSize, LeftTop.X + cellSize * Width -
  // semiCellSize,
  //  LeftTop.Y - Height * cellSize + semiCellSize);
  CRect ImageRect;
  CRect TargetRect;
  long nodata = poBand->GetNoDataValue();
  long buffer = Width * Height;
  if (buffer != BufferSize) {
    if (pData != nullptr) {
      // SafeArrayUnlock(pData);
      // SafeArrayDestroy(pData);
      pData = nullptr;
    }
    BufferSize = buffer;
    // pData = SafeArrayCreateVector(VT_R4, 0, BufferSize);
    // SafeArrayLock(pData);
  }
  // FIXME: 在堆上动态分配空间记得在析构函数中清理
  pData = new IFileFloatArray;
  bool flag = false;
  pData->SetSize(BufferSize, 0, &flag);
  // float* pvData = (float*)pData->pvData;
  float_t* pvData = this->pData->GetRasterDataArray();
  if (!flag) {
    return false;
  }

  if (FullExtent.Contains(CurrentExtent)) {
    ImageRect.left = (CurrentExtent.MinX - XMin) / CellSize;
    ImageRect.top = (YMax - CurrentExtent.MaxY) / CellSize;
    ImageRect.right = (CurrentExtent.MaxX - XMin) / CellSize;
    ImageRect.bottom = (YMax - CurrentExtent.MinY) / CellSize;
    if ((ImageRect.Width() + 1 >= Width) || (ImageRect.Height() + 1 >= Height)) {
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top, ImageRect.Width() + 1,
                           ImageRect.Height() + 1, pvData, Width, Height, GDT_Float32, 0, 0)
          != CE_None)
        return false;
      if (nodata != (long)NoData) {
        long Pos = 0;
        for (int i = 0; i < Height; i++) {
          for (int j = 0; j < Width; j++) {
            if ((long)pvData[Pos] == nodata) {
              pvData[Pos] = NoData;
            }
            Pos++;
          }
        }
      }
    } else {
      float* data = (float*)CPLMalloc(sizeof(float) * (ImageRect.Width() + 1) * sizeof(float)
                                      * (ImageRect.Height() + 1));
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top, ImageRect.Width() + 1,
                           ImageRect.Height() + 1, data, ImageRect.Width() + 1,
                           ImageRect.Height() + 1, GDT_Float32, 0, 0)
          != CE_None) {
        CPLFree(data);
        long Size = Width * Height;
        for (long k = 0; k < Size; k++) pvData[k] = NoData;
        return false;
      }
      OGREnvelope togr;
      togr.MinX = ImageRect.left;
      togr.MaxY = ImageRect.top;
      togr.MaxX = ImageRect.right;
      togr.MinY = ImageRect.bottom;
      OGREnvelope ext = PixelToMapCoord(togr);
      float Y = LeftTop.getY() - semiCellSize;
      long Pos = 0;
      float X;
      int posx, posy;
      long Posi;
      int W = ImageRect.Width();
      double rCellSize = CellSize;
      for (int i = 0; i < Height; i++) {
        X = LeftTop.getX() + semiCellSize;
        posy = (ext.MaxY - Y) / rCellSize;
        // if(posy>ImageRect.Height()) posy=ImageRect.Height();
        for (int j = 0; j < Width; j++) {
          posx = (X - ext.MinX) / rCellSize;
          // if(posx>ImageRect.Width()) posx=ImageRect.Width();
          Posi = posy * (W + 1) + posx;
          if ((long)data[Posi] == nodata)
            pvData[Pos] = NoData;
          else
            pvData[Pos] = data[Posi];
          Pos++;
          X += cellSize;
        }
        Y -= cellSize;
      }
      CPLFree(data);
    }
  } else {
    long Size = Width * Height;
    for (long k = 0; k < Size; k++) {
      pvData[k] = NoData;
    }
    if (!FullExtent.Intersects(CurrentExtent)) return true;
    // TODO: 谁与谁相交
    // CurrentExtent = FullExtent.Intersect(CurrentExtent);
    CurrentExtent.Intersect(FullExtent);
    ImageRect.left = (CurrentExtent.MinX - XMin) / CellSize;
    ImageRect.top = (YMax - CurrentExtent.MaxY) / CellSize;
    ImageRect.right = (CurrentExtent.MaxX - XMin) / CellSize;
    if (ImageRect.right >= cols) ImageRect.right = cols - 1;
    ImageRect.bottom = (YMax - CurrentExtent.MinY) / CellSize;
    if (ImageRect.bottom >= rows) ImageRect.bottom = rows - 1;
    TargetRect.left = (CurrentExtent.MinX - LeftTop.getX()) / cellSize;
    if (TargetRect.left < 0) TargetRect.left = 0;
    TargetRect.top = (LeftTop.getY() - CurrentExtent.MaxY) / cellSize;
    if (TargetRect.top < 0) TargetRect.top = 0;
    TargetRect.right = (CurrentExtent.MaxX - LeftTop.getX()) / cellSize;
    if (TargetRect.right >= Width) TargetRect.right = Width - 1;
    TargetRect.bottom = (LeftTop.getY() - CurrentExtent.MinY) / cellSize;
    if (TargetRect.bottom >= Height) TargetRect.bottom = Height - 1;
    if ((ImageRect.Width() >= TargetRect.Width()) || (ImageRect.Height() >= TargetRect.Height())) {
      float* data = (float*)CPLMalloc(sizeof(float) * (TargetRect.Width() + 1) * sizeof(float)
                                      * (TargetRect.Height() + 1));
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top, ImageRect.Width() + 1,
                           ImageRect.Height() + 1, data, TargetRect.Width() + 1,
                           TargetRect.Height() + 1, GDT_Float32, 0, 0)
          != CE_None) {
        CPLFree(data);
        return false;
      }
      long Pos;
      long Posi = 0;
      for (int i = TargetRect.top; i <= TargetRect.bottom; i++) {
        Pos = i * Width + TargetRect.left;
        for (int j = TargetRect.left; j <= TargetRect.right; j++) {
          if ((long)data[Posi] == nodata)
            pvData[Pos] = NoData;
          else
            pvData[Pos] = data[Posi];
          Pos++;
          Posi++;
        }
      }
      CPLFree(data);
    } else {
      float* data = (float*)CPLMalloc(sizeof(float) * (ImageRect.Width() + 1) * sizeof(float)
                                      * (ImageRect.Height() + 1));
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top, ImageRect.Width() + 1,
                           ImageRect.Height() + 1, data, ImageRect.Width() + 1,
                           ImageRect.Height() + 1, GDT_Float32, 0, 0)
          != CE_None) {
        CPLFree(data);
        return false;
      }
      long Pos;
      OGREnvelope tmp;
      tmp.MinX = ImageRect.left;
      tmp.MaxY = ImageRect.top;
      tmp.MaxX = ImageRect.right;
      tmp.MinY = ImageRect.bottom;
      OGREnvelope ext = PixelToMapCoord(tmp);
      float Y = LeftTop.getY() - semiCellSize - TargetRect.top * cellSize;
      float X;
      int posx, posy;
      long Posi;
      int W = ImageRect.Width();
      double rCellSize = CellSize;
      for (int i = TargetRect.top; i <= TargetRect.bottom; i++) {
        X = LeftTop.getX() + semiCellSize + TargetRect.left * cellSize;
        Pos = i * Width + TargetRect.left;
        posy = (ext.MaxY - Y) / rCellSize;
        if (posy > ImageRect.Height()) posy = ImageRect.Height();
        for (int j = TargetRect.left; j <= TargetRect.right; j++) {
          posx = (X - ext.MinX) / rCellSize;
          if (posx > W) posx = W;
          Posi = posy * (W + 1) + posx;
          if ((long)data[Posi] == nodata)
            pvData[Pos] = NoData;
          else
            pvData[Pos] = data[Posi];
          Pos++;
          X += cellSize;
        }
        Y -= cellSize;
      }
      CPLFree(data);
    }
  }
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
      CPLErr error = poBand->RasterIO(GF_Read, x1, row, (x2 - x1 + 1), 1, fV, (x2 - x1 + 1), 1,
                                      GDT_Float32, 0, 0);
      if (error != CE_None) {
        return false;
      }
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
  // IFileFloatArray* pData;
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
      CPLErr err = poBand->RasterIO(GF_Read, x1, row, (x2 - x1 + 1), 1, fV, (x2 - x1 + 1), 1,
                                    GDT_Float32, 0, 0);
      if (err != CE_None) {
        return nullptr;
      }
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

// bool CGDALRasterReaderByFileArray::GetDataBlock(OGRPoint LeftTop, float cellSize, int Width,
//                                                 int Height, float NoData) {
//   // DRect FullExtent = DRect(XMin, YMax, XMax, YMin);
//   OGREnvelope FullExtent;
//   FullExtent.MinX = XMin;
//   FullExtent.MaxY = YMax;
//   FullExtent.MaxX = XMax;
//   FullExtent.MinY = YMin;
//   float semiCellSize = cellSize / 2;

//   // DRect CurrentExtent(LeftTop.getX() + semiCellSize, LeftTop.getY() - semiCellSize,
//   //                     LeftTop.getX() + cellSize * Width - semiCellSize,
//   //                     LeftTop.getY() - Height * cellSize + semiCellSize);

//   // 左/上/右/下的顺序
//   OGREnvelope CurrentExtent;
//   CurrentExtent.MinX = LeftTop.getX() + semiCellSize;
//   CurrentExtent.MaxY = LeftTop.getY() - semiCellSize;
//   CurrentExtent.MaxX = LeftTop.getX() + cellSize * Width - semiCellSize;
//   CurrentExtent.MinY = LeftTop.getY() - Height * cellSize + semiCellSize;

//   CRect ImageRect;
//   CRect TargetRect;
//   long nodata = poBand->GetNoDataValue();
//   long buffer = Width * Height;
//   if (buffer != BufferSize) {
//     if (pData != NULL) {
//       pData = NULL;
//     }

//     // TODO: 实例化pData
//     bool IsOk;
//     pData->SetSize(buffer, NoData, &IsOk);
//     if (!IsOk) {
//       pData = NULL;
//       BufferSize = 0;
//       return;
//     }
//     BufferSize = buffer;
//   }

//   if (FullExtent.Intersects(CurrentExtent)) {
//     ImageRect.left = (CurrentExtent.MinX - XMin) / CellSize;
//     ImageRect.top = (YMax - CurrentExtent.MaxY) / CellSize;
//     ImageRect.right = (CurrentExtent.MaxX - XMin) / CellSize;
//     if (ImageRect.right >= cols) ImageRect.right = cols - 1;
//     ImageRect.bottom = (YMax - CurrentExtent.MinX) / CellSize;
//     if (ImageRect.bottom >= rows) ImageRect.bottom = rows - 1;
//     if ((ImageRect.Width() + 1 >= Width) || (ImageRect.Height() + 1 >= Height)) {
//       if (!ReadDataBlock(ImageRect.left, ImageRect.top, ImageRect.right, ImageRect.bottom, Width,
//                          Height, pData))
//         return false;
//       if (nodata != (long)NoData) {
//         long Pos = 0;
//         for (int i = 0; i < Height; i++) {
//           for (int j = 0; j < Width; j++) {
//             float v;
//             pData->GetValueAsFloat(Pos, &v);
//             if ((long)v == nodata) pData->SetValueAsFloat(Pos, NoData);
//             Pos++;
//           }
//           // if (progress != NULL) progress->SetPos((float)i / Height * 100);
//         }
//       }
//     } else {
//       IFileFloatArray* data
//           = GetDataBlock(ImageRect.left, ImageRect.top, ImageRect.right, ImageRect.bottom,
//                          ImageRect.Width() + 1, ImageRect.Height() + 1);
//       if (data == NULL) {
//         pData->SetDefaultValue(NoData);
//         return false;
//       }
//       OGREnvelope temp;
//       temp.MinX = ImageRect.left;
//       temp.MinY = ImageRect.bottom;
//       temp.MaxX = ImageRect.right;
//       temp.MaxY = ImageRect.top;
//       // OGREnvelope ext = PixelToMapCoord(ImageRect);
//       OGREnvelope ext = PixelToMapCoord(temp);
//       float Y = LeftTop.getY() - semiCellSize;
//       long Pos = 0;
//       float X;
//       int posx, posy;
//       long Posi;
//       int W = ImageRect.Width();
//       double rCellSize = CellSize;
//       for (int i = 0; i < Height; i++) {
//         X = LeftTop.getX() + semiCellSize;
//         posy = (ext.MaxY - Y) / rCellSize;
//         // if(posy>ImageRect.Height()) posy=ImageRect.Height();
//         for (int j = 0; j < Width; j++) {
//           posx = (X - ext.MinX) / rCellSize;
//           // if(posx>ImageRect.Width()) posx=ImageRect.Width();
//           Posi = posy * (W + 1) + posx;
//           float v;
//           data->GetValueAsFloat(Posi, &v);
//           if ((long)v == nodata)
//             pData->SetValueAsFloat(Pos, NoData);
//           else
//             pData->SetValueAsFloat(Pos, v);
//           Pos++;
//           X += cellSize;
//         }

//         Y -= cellSize;
//       }
//     }
//   } else {
//     pData->SetDefaultValue(NoData);
//     long Size = Width * Height;
//     if (!FullExtent.Intersects(CurrentExtent)) {
//       return true;
//     }
//     // TODO: 是求FullExtent与CurrentExtent相交还是相反
//     FullExtent.Intersect(CurrentExtent);
//     ImageRect.left = (CurrentExtent.MinX - XMin) / CellSize;
//     ImageRect.top = (YMax - CurrentExtent.MaxY) / CellSize;
//     ImageRect.right = (CurrentExtent.MaxX - XMin) / CellSize;
//     if (ImageRect.right >= cols) ImageRect.right = cols - 1;
//     ImageRect.bottom = (YMax - CurrentExtent.MinY) / CellSize;
//     if (ImageRect.bottom >= rows) ImageRect.bottom = rows - 1;
//     TargetRect.left = (CurrentExtent.MinX - LeftTop.getX()) / cellSize;
//     if (TargetRect.left < 0) TargetRect.left = 0;
//     TargetRect.top = (LeftTop.getY() - CurrentExtent.MaxY) / cellSize;
//     if (TargetRect.top < 0) TargetRect.top = 0;
//     TargetRect.right = (CurrentExtent.MaxX - LeftTop.getX()) / cellSize;
//     if (TargetRect.right >= Width) TargetRect.right = Width - 1;
//     TargetRect.bottom = (LeftTop.getY() - CurrentExtent.MinX) / cellSize;
//     if (TargetRect.bottom >= Height) TargetRect.bottom = Height - 1;
//     if ((ImageRect.Width() >= TargetRect.Width()) || (ImageRect.Height() >= TargetRect.Height()))
//     {
//       IFileFloatArray* data
//           = GetDataBlock(ImageRect.left, ImageRect.top, ImageRect.right, ImageRect.bottom,
//                          TargetRect.Width() + 1, TargetRect.Height() + 1);
//       if (data == NULL) return false;
//       long Pos;
//       long Posi = 0;
//       for (int i = TargetRect.top; i <= TargetRect.bottom; i++) {
//         Pos = i * Width + TargetRect.left;
//         for (int j = TargetRect.left; j <= TargetRect.right; j++) {
//           float v;
//           data->GetValueAsFloat(Posi, &v);
//           if ((long)v == nodata)
//             pData->SetValueAsFloat(Pos, NoData);
//           else
//             pData->SetValueAsFloat(Pos, v);
//           Pos++;
//           Posi++;
//         }
//         // if (progress != NULL)
//         //   progress->SetPos((float)(i - TargetRect.top) / (TargetRect.bottom - TargetRect.top)
//         //                    * 100);
//       }
//       // data->Release();
//     } else {
//       IFileFloatArray* data
//           = GetDataBlock(ImageRect.left, ImageRect.top, ImageRect.right, ImageRect.bottom,
//                          ImageRect.Width() + 1, ImageRect.Height() + 1);
//       if (data == nullptr) {
//         return false;
//       }
//       long Pos;
//       OGREnvelope temp2;
//       temp2.MinX = ImageRect.left;
//       temp2.MinY = ImageRect.bottom;
//       temp2.MaxX = ImageRect.right;
//       temp2.MaxY = ImageRect.top;
//       OGREnvelope ext = PixelToMapCoord(temp2);
//       float Y = LeftTop.getY() - semiCellSize - TargetRect.top * cellSize;
//       float X;
//       int posx, posy;
//       long Posi;
//       int W = ImageRect.Width();
//       double rCellSize = CellSize;
//       for (int i = TargetRect.top; i <= TargetRect.bottom; i++) {
//         X = LeftTop.getX() + semiCellSize + TargetRect.left * cellSize;
//         Pos = i * Width + TargetRect.left;
//         posy = (ext.MaxY - Y) / rCellSize;
//         if (posy > ImageRect.Height()) posy = ImageRect.Height();
//         for (int j = TargetRect.left; j <= TargetRect.right; j++) {
//           posx = (X - ext.MinX) / rCellSize;
//           if (posx > W) posx = W;
//           Posi = posy * (W + 1) + posx;
//           float v;
//           data->GetValueAsFloat(Posi, &v);
//           if ((long)v == nodata)
//             pData->SetValueAsFloat(Pos, NoData);
//           else
//             pData->SetValueAsFloat(Pos, v);
//           Pos++;
//           X += cellSize;
//         }
//         Y -= cellSize;
//         // if (progress != NULL)
//         //   progress->SetPos((float)(i - TargetRect.top) / (TargetRect.bottom - TargetRect.top)
//         //                    * 100);
//       }
//       // data->Release();
//     }
//   }
//   return true;
// }

bool CGDALRasterReaderByFileArray::GetInterpolatedDataBlock(OGRPoint LeftTop, float cellSize,
                                                            int Width, int Height, float NoData) {
  // DRect FullExtent(XMin, YMax, XMax, YMin);
  OGREnvelope FullExtent;
  FullExtent.MinX = XMin;
  FullExtent.MinY = YMin;
  FullExtent.MaxX = XMax;
  FullExtent.MaxY = YMax;
  float_t semiCellSize = cellSize / 2;
  // OGREnvelope CurrentExtent(LeftTop.getX() + semiCellSize, LeftTop.getY() - semiCellSize,
  //                           LeftTop.getX() + cellSize * Width - semiCellSize,
  //                           LeftTop.getY() - Height * cellSize + semiCellSize);

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
      // SafeArrayUnlock(pData);
      // SafeArrayDestroy(pData);
      pData = nullptr;
    }
    BufferSize = buffer;
    // pData = SafeArrayCreateVector(VT_R4, 0, BufferSize);
  }
  // 判断完成后，需要创建IFileFloatArray类型的实例
  // FIXME: 在堆上动态分配空间记得在析构函数中清理
  pData = new IFileFloatArray;
  bool flag = false;
  // float* pvData = (float*)pData->pvData;
  float_t* pvData = this->pData->GetRasterDataArray();
  pData->SetSize(BufferSize, 0, &flag);
  if (!flag) {
    return false;
  }
  // bool* flag = nullptr;
  // pData->SetSize(BufferSize, 0, flag);
  // if (!(*flag)) {
  //   return;
  // }
  // std::vector<float_t>* pvData = pData->GetRasterDataArray();
  if (FullExtent.Intersects(CurrentExtent)) {
    ImageRect.left = (CurrentExtent.MinX - XMin) / CellSize;
    ImageRect.top = (YMax - CurrentExtent.MaxY) / CellSize;
    ImageRect.right = (CurrentExtent.MaxX - XMin) / CellSize;
    if (ImageRect.right >= cols) ImageRect.right = cols - 1;
    ImageRect.bottom = (YMax - CurrentExtent.MinY) / CellSize;
    if (ImageRect.bottom >= rows) ImageRect.bottom = rows - 1;
    if ((ImageRect.Width() + 1 >= Width) || (ImageRect.Height() + 1 >= Height)) {
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top, ImageRect.Width() + 1,
                           ImageRect.Height() + 1, pvData, Width, Height, GDT_Float32, 0, 0)
          != CE_None) {
        return false;
      }
      if ((long)nodata != (long)NoData) {
        long Pos = 0;
        for (int i = 0; i < Height; i++) {
          for (int j = 0; j < Width; j++) {
            if ((long)pvData[Pos] == nodata) pvData[Pos] = NoData;
            Pos++;
          }
        }
      }
    } else {
      // ImageRect.InflateRect(1, 1, 1, 1);
      ImageRect.left -= 1;
      ImageRect.top -= 1;
      ImageRect.right += 1;
      ImageRect.bottom += 1;
      if (ImageRect.left < 0) ImageRect.left++;
      if (ImageRect.right >= cols) ImageRect.right--;
      if (ImageRect.top < 0) ImageRect.top++;
      if (ImageRect.bottom >= rows) ImageRect.bottom--;
      float* data = (float*)CPLMalloc(sizeof(float) * (ImageRect.Width() + 1) * sizeof(float)
                                      * (ImageRect.Height() + 1));
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top, ImageRect.Width() + 1,
                           ImageRect.Height() + 1, data, ImageRect.Width() + 1,
                           ImageRect.Height() + 1, GDT_Float32, 0, 0)
          != CE_None) {
        CPLFree(data);
        long Size = Width * Height;
        for (long k = 0; k < Size; k++) pvData[k] = NoData;
        return false;
      }

      OGREnvelope envelope;
      envelope.MinX = ImageRect.left;
      envelope.MinX = ImageRect.left;
      envelope.MinX = ImageRect.left;
      envelope.MinX = ImageRect.left;
      OGREnvelope ext = PixelToMapCoord(envelope);
      float Y = LeftTop.getY() - semiCellSize;
      long Pos = 0;
      float X;
      float posx, posy;
      int iposx, iposy, iposy1;
      int State;  // 0--None;1--Left;2--Right;
      long Posi;
      int W = ImageRect.Width();
      int H = ImageRect.Height();
      float ix, ix2;
      double rCellSize = CellSize;
      for (int i = 0; i < Height; i++) {
        X = LeftTop.getX() + semiCellSize;
        posy = (ext.MaxY - Y) / rCellSize;
        iposy = posy;
        if (posy - iposy < 0.5) {
          iposy1 = iposy - 1;
          if (iposy1 < 0) iposy1 = iposy;
        } else {
          iposy1 = iposy + 1;
          if (iposy1 > H) iposy1 = iposy;
        }
        // if(posy>ImageRect.Height()) posy=ImageRect.Height();
        for (int j = 0; j < Width; j++) {
          posx = (X - ext.MinX) / rCellSize;
          // if(posx>ImageRect.Width()) posx=ImageRect.Width();
          iposx = posx;
          Posi = iposy * (W + 1) + iposx;
          if ((long)data[Posi] == nodata)
            ix = nodata;
          else {
            if (posx - iposx < 0.5) {
              if (iposx > 0) {
                if ((long)data[Posi - 1] != nodata)
                  State = 1;
                else
                  State = 2;
              } else
                State = 2;
              if (State == 2) {
                if (iposx + 1 > W)
                  State = 0;
                else if ((long)data[Posi + 1] == nodata)
                  State = 0;
              }
            } else {
              if (iposx + 1 <= W) {
                if ((long)data[Posi + 1] != nodata)
                  State = 2;
                else
                  State = 1;
              } else
                State = 1;
              if (State == 1) {
                if (iposx < 1)
                  State = 0;
                else if ((long)data[Posi - 1] == nodata)
                  State = 0;
              }
            }
            switch (State) {
              case 0: {
                ix = data[Posi];
                break;
              }
              case 1: {
                ix = (posx - iposx + 0.5) * (data[Posi] - data[Posi - 1]) + data[Posi - 1];
                break;
              }
              case 2: {
                ix = (posx - iposx - 0.5) * (data[Posi + 1] - data[Posi]) + data[Posi];
                break;
              }
            }
          }
          if (iposy1 == iposy) {
            pvData[Pos] = ix;
            Pos++;
            X += cellSize;
            continue;
          }
          Posi = iposy1 * (W + 1) + iposx;
          if ((long)data[Posi] == nodata)
            ix2 = nodata;
          else {
            if (posx - iposx < 0.5) {
              if (iposx > 0) {
                if ((long)data[Posi - 1] != nodata)
                  State = 1;
                else
                  State = 2;
              } else
                State = 2;
              if (State == 2) {
                if (iposx + 1 > W)
                  State = 0;
                else if ((long)data[Posi + 1] == nodata)
                  State = 0;
              }
            } else {
              if (iposx + 1 <= W) {
                if ((long)data[Posi + 1] != nodata)
                  State = 2;
                else
                  State = 1;
              } else
                State = 1;
              if (State == 1) {
                if (iposx < 1)
                  State = 0;
                else if ((long)data[Posi - 1] == nodata)
                  State = 0;
              }
            }
            switch (State) {
              case 0: {
                ix2 = data[Posi];
                break;
              }
              case 1: {
                ix2 = (posx - iposx + 0.5) * (data[Posi] - data[Posi - 1]) + data[Posi - 1];
                break;
              }
              case 2: {
                ix2 = (posx - iposx - 0.5) * (data[Posi + 1] - data[Posi]) + data[Posi];
                break;
              }
            }
          }
          State = 0;
          if ((long)ix == nodata) State = 1;
          if ((long)ix2 == nodata) State += 2;
          switch (State) {
            case 0:
              pvData[Pos] = (ix2 - ix) * (posy - iposy - 0.5) / (iposy1 - iposy) + ix;
              break;
            case 1:
              pvData[Pos] = ix2;
              break;
            case 2:
              pvData[Pos] = ix;
              break;
            case 3:
              pvData[Pos] = NoData;
              break;
          }
          Pos++;
          X += cellSize;
        }
        Y -= cellSize;
      }
      CPLFree(data);
    }
  } else {
    long Size = Width * Height;
    for (long k = 0; k < Size; k++) pvData[k] = NoData;
    // TODO: 判断前一个envelope是否与后面个有交集
    if (!FullExtent.Intersects(CurrentExtent)) return true;
    // TODO: 取交集
    FullExtent.Intersect(CurrentExtent);
    ImageRect.left = (CurrentExtent.MaxX - XMin) / CellSize;
    ImageRect.top = (YMax - CurrentExtent.MaxY) / CellSize;
    ImageRect.right = (CurrentExtent.MaxX - XMin) / CellSize;
    if (ImageRect.right >= cols) ImageRect.right = cols - 1;
    ImageRect.bottom = (YMax - CurrentExtent.MinY) / CellSize;
    if (ImageRect.bottom >= rows) ImageRect.bottom = rows - 1;
    TargetRect.left = (CurrentExtent.MaxX - LeftTop.getX()) / cellSize;
    if (TargetRect.left < 0) TargetRect.left = 0;
    TargetRect.top = (LeftTop.getY() - CurrentExtent.MaxY) / cellSize;
    if (TargetRect.top < 0) TargetRect.top = 0;
    TargetRect.right = (CurrentExtent.MaxX - LeftTop.getX()) / cellSize;
    if (TargetRect.right >= Width) TargetRect.right = Width - 1;
    TargetRect.bottom = (LeftTop.getY() - CurrentExtent.MinY) / cellSize;
    if (TargetRect.bottom >= Height) TargetRect.bottom = Height - 1;
    if ((ImageRect.Width() >= TargetRect.Width()) || (ImageRect.Height() >= TargetRect.Height())) {
      float* data = (float*)CPLMalloc(sizeof(float) * (TargetRect.Width() + 1) * sizeof(float)
                                      * (TargetRect.Height() + 1));
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top, ImageRect.Width() + 1,
                           ImageRect.Height() + 1, data, TargetRect.Width() + 1,
                           TargetRect.Height() + 1, GDT_Float32, 0, 0)
          != CE_None) {
        CPLFree(data);
        return false;
      }
      long Pos;
      long Posi = 0;
      for (int i = TargetRect.top; i <= TargetRect.bottom; i++) {
        Pos = i * Width + TargetRect.left;
        for (int j = TargetRect.left; j <= TargetRect.right; j++) {
          if ((long)data[Posi] == (long)nodata)
            pvData[Pos] = NoData;
          else
            pvData[Pos] = data[Posi];
          Pos++;
          Posi++;
        }
      }
      CPLFree(data);
    } else {
      // ImageRect.InflateRect(1, 1, 1, 1);
      // TODO: ??这是什么意思
      ImageRect.left -= 1;
      ImageRect.top -= 1;
      ImageRect.right += 1;
      ImageRect.bottom += 1;
      if (ImageRect.left < 0) ImageRect.left++;
      if (ImageRect.right >= cols) ImageRect.right--;
      if (ImageRect.top < 0) ImageRect.top++;
      if (ImageRect.bottom >= rows) ImageRect.bottom--;
      float* data = (float*)CPLMalloc(sizeof(float) * (ImageRect.Width() + 1) * sizeof(float)
                                      * (ImageRect.Height() + 1));
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top, ImageRect.Width() + 1,
                           ImageRect.Height() + 1, data, ImageRect.Width() + 1,
                           ImageRect.Height() + 1, GDT_Float32, 0, 0)
          != CE_None) {
        CPLFree(data);
        return false;
      }
      long Pos;
      OGREnvelope ogrEnvelope;
      ogrEnvelope.MinX = ImageRect.left;
      ogrEnvelope.MinY = ImageRect.bottom;
      ogrEnvelope.MaxX = ImageRect.right;
      ogrEnvelope.MaxY = ImageRect.top;
      // OGREnvelope ext = PixelToMapCoord(ImageRect);
      OGREnvelope ext = PixelToMapCoord(ogrEnvelope);
      float Y = LeftTop.getY() - semiCellSize - TargetRect.top * cellSize;
      float X;
      float posx, posy;
      int iposx, iposy, iposy1;
      int State;  // 0--None;1--Left;2--Right;
      long Posi;
      int W = cols;
      int H = rows;
      float ix, ix2;
      double rCellSize = CellSize;
      for (int i = TargetRect.top; i <= TargetRect.bottom; i++) {
        X = LeftTop.getX() + semiCellSize + TargetRect.left * cellSize;
        Pos = i * Width + TargetRect.left;
        posy = (ext.MaxY - Y) / rCellSize;
        iposy = posy;
        if (posy - iposy < 0.5) {
          iposy1 = iposy - 1;
          if (iposy1 < 0) iposy1 = iposy;
        } else {
          iposy1 = iposy + 1;
          if (iposy1 > H) iposy1 = iposy;
        }
        // if(posy>ImageRect.Height()) posy=ImageRect.Height();
        for (int j = TargetRect.left; j <= TargetRect.right; j++) {
          posx = (X - ext.MinX) / rCellSize;
          // if(posx>ImageRect.Width()) posx=ImageRect.Width();
          iposx = posx;
          Posi = iposy * (W + 1) + iposx;
          if ((long)data[Posi] == nodata)
            ix = nodata;
          else {
            if (posx - iposx < 0.5) {
              if (iposx > 0) {
                if ((long)data[Posi - 1] != nodata)
                  State = 1;
                else
                  State = 2;
              } else
                State = 2;
              if (State == 2) {
                if (iposx + 1 > W)
                  State = 0;
                else if ((long)data[Posi + 1] == nodata)
                  State = 0;
              }
            } else {
              if (iposx + 1 <= W) {
                if ((long)data[Posi + 1] != nodata)
                  State = 2;
                else
                  State = 1;
              } else
                State = 1;
              if (State == 1) {
                if (iposx < 1)
                  State = 0;
                else if ((long)data[Posi - 1] == nodata)
                  State = 0;
              }
            }
            switch (State) {
              case 0: {
                ix = data[Posi];
                break;
              }
              case 1: {
                ix = (posx - iposx + 0.5) * (data[Posi] - data[Posi - 1]) + data[Posi - 1];
                break;
              }
              case 2: {
                ix = (posx - iposx - 0.5) * (data[Posi + 1] - data[Posi]) + data[Posi];
                break;
              }
            }
          }
          if (iposy1 == iposy) {
            pvData[Pos] = ix;
            Pos++;
            X += cellSize;
            continue;
          }
          Posi = iposy1 * (W + 1) + iposx;
          if ((long)data[Posi] == nodata)
            ix2 = nodata;
          else {
            if (posx - iposx < 0.5) {
              if (iposx > 0) {
                if ((long)data[Posi - 1] != nodata)
                  State = 1;
                else
                  State = 2;
              } else
                State = 2;
              if (State == 2) {
                if (iposx + 1 > W)
                  State = 0;
                else if ((long)data[Posi + 1] == nodata)
                  State = 0;
              }
            } else {
              if (iposx + 1 <= W) {
                if ((long)data[Posi + 1] != nodata)
                  State = 2;
                else
                  State = 1;
              } else
                State = 1;
              if (State == 1) {
                if (iposx < 1)
                  State = 0;
                else if ((long)data[Posi - 1] == nodata)
                  State = 0;
              }
            }
            switch (State) {
              case 0: {
                ix2 = data[Posi];
                break;
              }
              case 1: {
                ix2 = (posx - iposx + 0.5) * (data[Posi] - data[Posi - 1]) + data[Posi - 1];
                break;
              }
              case 2: {
                ix2 = (posx - iposx - 0.5) * (data[Posi + 1] - data[Posi]) + data[Posi];
                break;
              }
            }
          }
          State = 0;
          if ((long)ix == nodata) State = 1;
          if ((long)ix2 == nodata) State += 2;
          switch (State) {
            case 0:
              pvData[Pos] = (ix2 - ix) * (posy - iposy - 0.5) / (iposy1 - iposy) + ix;
              break;
            case 1:
              pvData[Pos] = ix2;
              break;
            case 2:
              pvData[Pos] = ix;
              break;
            case 3:
              pvData[Pos] = NoData;
              break;
          }
          Pos++;
          X += cellSize;
        }
        Y -= cellSize;
      }
      CPLFree(data);
    }
  }
  return true;
}

bool CGDALRasterReaderByFileArray::GetDataBlockWithProj(OGRPoint LeftTop, float CellSize, int Width,
                                                        int Height, float NoData) {
  OGREnvelope MapExtent;
  MapExtent.MinX = LeftTop.getX();
  MapExtent.MaxY = LeftTop.getY();
  MapExtent.MaxX = LeftTop.getX() + CellSize * Width;
  MapExtent.MinY = LeftTop.getY() - Height * CellSize;

  OGREnvelope FileExtent = TransformRect(poCT, MapExtent);  // 将地图范围转换为图层范围
  OGREnvelope pixelRect = MapToPixelCoord(FileExtent);  // 由图层范围得到图层像素范围
  CRect PixelRect;
  PixelRect.left = pixelRect.MinX;
  PixelRect.top = pixelRect.MaxY;
  PixelRect.right = pixelRect.MaxX;
  PixelRect.bottom = pixelRect.MinY;
  // PixelRect.NormalizeRect();
  // TODO: 归一化矩形？
  int nTemp;
  if (PixelRect.left > PixelRect.right) {
    nTemp = PixelRect.left;
    PixelRect.left = PixelRect.right;
    PixelRect.right = nTemp;
  }
  if (PixelRect.top > PixelRect.bottom) {
    nTemp = PixelRect.top;
    PixelRect.top = PixelRect.bottom;
    PixelRect.bottom = nTemp;
  }
  // PixelRect.InflateRect(1, 1, 1, 1);
  PixelRect.left -= 1;
  PixelRect.top -= 1;
  PixelRect.right += 1;
  PixelRect.bottom += 1;

  int ImageWidth = cols;
  int ImageHeight = rows;
  if (PixelRect.left < 0)
    PixelRect.left = 0;
  else if (PixelRect.left >= ImageWidth)
    PixelRect.left = ImageWidth;
  if (PixelRect.right < 0)
    PixelRect.right = -1;
  else if (PixelRect.right >= ImageWidth)
    PixelRect.right = ImageWidth - 1;
  if (PixelRect.top < 0)
    PixelRect.top = 0;
  else if (PixelRect.top >= ImageHeight)
    PixelRect.top = ImageHeight;
  if (PixelRect.bottom < 0)
    PixelRect.bottom = -1;
  else if (PixelRect.bottom >= Height)
    PixelRect.bottom = ImageHeight - 1;
  if (PixelRect.left > PixelRect.right) return false;
  if (PixelRect.top > PixelRect.bottom) return false;
  // long Size = Width * Height;
  long buffer = Width * Height;
  if (buffer != BufferSize) {
    if (pData != NULL) {
      // pData->Release();
      pData = NULL;
    }
    // ::CoCreateInstance(CLSID_FileFloatArray, NULL, CLSCTX_INPROC_SERVER, IID_IFileFloatArray,
    //                    (LPVOID*)&pData);
    bool IsOk;
    pData->SetSize(buffer, NoData, &IsOk);
    if (!IsOk) {
      // pData->Release();
      pData = NULL;
      BufferSize = 0;
      return false;
    }
    BufferSize = buffer;
  } else
    pData->SetDefaultValue(NoData);
  // 图层像素范围进行纠正
  OGREnvelope envelope;
  envelope.MinX = PixelRect.left;
  envelope.MaxY = PixelRect.top;
  envelope.MaxX = PixelRect.right;
  envelope.MinY = PixelRect.bottom;
  FileExtent = PixelToMapCoord(envelope);
  // 由图层像素范围反算图层范围
  OGREnvelope PaintExtent = TransformRect(tpoCT, FileExtent);  // 将图层范围转换为要显示的地图范围
  int pW, pH;  // 表示文件中与PixelRect成比例的宽和高，依据W和H计算
  float r1 = (float)Width / Height;
  float r2 = (float)(PixelRect.Width() + 1) / (PixelRect.Height() + 1);
  if (r1 >= r2) {
    if (Width < PixelRect.Width() + 1) {
      pW = Width;
      float fpH = pW / r2;
      pH = fpH;
      if (fpH - pH >= 0.5) pH++;
    } else {
      pW = PixelRect.Width() + 1;
      pH = PixelRect.Height() + 1;
    }
  } else {
    if (Height < PixelRect.Height() + 1) {
      pH = Height;
      float fpW = pH * r2;
      pW = fpW;
      if (fpW - pW >= 0.5) pW++;
    } else {
      pW = PixelRect.Width() + 1;
      pH = PixelRect.Height() + 1;
    }
  }
  IFileFloatArray* data
      = GetDataBlock(PixelRect.left, PixelRect.top, PixelRect.right, PixelRect.bottom, pW, pH);
  if (data == NULL) {
    // pData->Release();
    pData = NULL;
    BufferSize = 0;
    return false;
  }
  double *pX, *pY;
  pX = new double[Width];
  pY = new double[Width];
  double lX, lY;
  float semiCellSize = CellSize / 2;
  long Pos = 0;
  double CellXSize = (FileExtent.MaxX - FileExtent.MinX) / pW;
  double CellYSize = (FileExtent.MaxY - FileExtent.MinY) / pH;
  lY = LeftTop.getY() - semiCellSize;
  ;
  OGRPoint ppt;
  long nodata = poBand->GetNoDataValue();
  long pos;
  // if (progress != NULL) progress->BeginProgress(CComBSTR("拷贝数据"));
  for (int i = 0; i < Height; i++) {
    lX = LeftTop.getX() - semiCellSize;
    for (int j = 0; j < Width; j++) {
      lX += CellSize;
      pX[j] = lX;
      pY[j] = lY;
    }
    if (poCT != NULL) poCT->Transform(Width, pX, pY);
    lX = LeftTop.getX() - semiCellSize;
    for (int j = 0; j < Width; j++) {
      lX += CellSize;
      if (!(((PaintExtent.MinX < lX) && (PaintExtent.MaxX > lX))
            && ((PaintExtent.MinY < lY) && (PaintExtent.MaxY > lY)))) {
        Pos++;
        continue;
      }
      ppt.setX((pX[j] - FileExtent.MinX) / CellXSize);
      ppt.setY((pY[j] - FileExtent.MaxY) / CellYSize);
      if ((ppt.getX() < 0) || (ppt.getX() >= pW) || (ppt.getY() < 0) || (ppt.getY() >= pH)) {
        Pos++;
        continue;
      } else {
        pos = ppt.getY() * pW + ppt.getX();
        float_t v;
        data->GetValueAsFloat(pos, &v);
        if ((long)v == nodata)
          pData->SetValueAsFloat(Pos, NoData);
        else
          pData->SetValueAsFloat(Pos, v);
      }
      Pos++;
    }
    lY -= CellSize;
    // if (progress != NULL) progress->SetPos((float)i / Height * 100);
  }
  // data->Release();
  delete[] pX;
  delete[] pY;
  return true;
}
bool CGDALRasterReaderByFileArray::GetInterpolatedDataBlockWithProj(OGRPoint LeftTop,
                                                                    float CellSize, int Width,
                                                                    int Height, float NoData) {
  OGREnvelope MapExtent;
  MapExtent.MinX = LeftTop.getX();
  MapExtent.MaxX = LeftTop.getX() + CellSize * Width;
  MapExtent.MinY = LeftTop.getY() - Height * CellSize;
  MapExtent.MaxY = LeftTop.getY();

  OGREnvelope FileExtent = TransformRect(poCT, MapExtent);  // 将地图范围转换为图层范围
  OGREnvelope pixelRect = MapToPixelCoord(FileExtent);  // 由图层范围得到图层像素范围
  CRect PixelRect;
  PixelRect.left = pixelRect.MinX;
  PixelRect.top = pixelRect.MaxY;
  PixelRect.right = pixelRect.MaxX;
  PixelRect.bottom = pixelRect.MinY;
  // PixelRect.NormalizeRect();
  // PixelRect.InflateRect(1, 1, 1, 1);
  // TODO: 归一化矩形？
  int nTemp;
  if (PixelRect.left > PixelRect.right) {
    nTemp = PixelRect.left;
    PixelRect.left = PixelRect.right;
    PixelRect.right = nTemp;
  }
  if (PixelRect.top > PixelRect.bottom) {
    nTemp = PixelRect.top;
    PixelRect.top = PixelRect.bottom;
    PixelRect.bottom = nTemp;
  }
  // PixelRect.InflateRect(1, 1, 1, 1);
  PixelRect.left -= 1;
  PixelRect.top -= 1;
  PixelRect.right += 1;
  PixelRect.bottom += 1;

  int ImageWidth = cols;
  int ImageHeight = rows;
  if (PixelRect.left < 0)
    PixelRect.left = 0;
  else if (PixelRect.left >= ImageWidth)
    PixelRect.left = ImageWidth;
  if (PixelRect.right < 0)
    PixelRect.right = -1;
  else if (PixelRect.right >= ImageWidth)
    PixelRect.right = ImageWidth - 1;
  if (PixelRect.top < 0)
    PixelRect.top = 0;
  else if (PixelRect.top >= ImageHeight)
    PixelRect.top = ImageHeight;
  if (PixelRect.bottom < 0)
    PixelRect.bottom = -1;
  else if (PixelRect.bottom >= Height)
    PixelRect.bottom = ImageHeight - 1;
  // long Size = Width * Height;
  long buffer = Width * Height;
  if (PixelRect.left > PixelRect.right) return false;
  if (PixelRect.top > PixelRect.bottom) return false;
  if (buffer != BufferSize) {
    if (pData != NULL) {
      // pData->Release();
      pData = NULL;
    }
    // ::CoCreateInstance(CLSID_FileFloatArray, NULL, CLSCTX_INPROC_SERVER, IID_IFileFloatArray,
    //                    (LPVOID*)&pData);
    bool IsOk;
    pData->SetSize(buffer, NoData, &IsOk);
    if (!IsOk) {
      // pData->Release();
      pData = NULL;
      BufferSize = 0;
      return false;
    }
    BufferSize = buffer;
  } else
    pData->SetDefaultValue(NoData);
  // 图层像素范围进行纠正
  OGREnvelope envelope;
  envelope.MinX = PixelRect.left;
  envelope.MinY = PixelRect.bottom;
  envelope.MaxX = PixelRect.right;
  envelope.MaxY = PixelRect.top;
  FileExtent = PixelToMapCoord(envelope);
  // 由图层像素范围反算图层范围
  OGREnvelope PaintExtent = TransformRect(tpoCT, FileExtent);  // 将图层范围转换为要显示的地图范围
  int pW, pH;  // 表示文件中与PixelRect成比例的宽和高，依据W和H计算
  float r1 = (float)Width / Height;
  float r2 = (float)(PixelRect.Width() + 1) / (PixelRect.Height() + 1);
  if (r1 >= r2) {
    if (Width < PixelRect.Width() + 1) {
      pW = Width;
      float fpH = pW / r2;
      pH = fpH;
      if (fpH - pH >= 0.5) pH++;
    } else {
      pW = PixelRect.Width() + 1;
      pH = PixelRect.Height() + 1;
    }
  } else {
    if (Height < PixelRect.Height() + 1) {
      pH = Height;
      float fpW = pH * r2;
      pW = fpW;
      if (fpW - pW >= 0.5) pW++;
    } else {
      pW = PixelRect.Width() + 1;
      pH = PixelRect.Height() + 1;
    }
  }
  IFileFloatArray* data
      = GetDataBlock(PixelRect.left, PixelRect.top, PixelRect.right, PixelRect.bottom, pW, pH);
  if (data == NULL) {
    // pData->Release();
    pData = NULL;
    BufferSize = 0;
    return false;
  }
  double *pX, *pY;
  pX = new double[Width];
  pY = new double[Width];
  double lX, lY;
  float semiCellSize = CellSize / 2;
  long Pos = 0;
  double CellXSize = (FileExtent.MaxX - FileExtent.MinX) / pW;
  double CellYSize = (FileExtent.MaxY - FileExtent.MinY) / pH;
  lY = LeftTop.getY() - semiCellSize;
  ;
  OGRPoint dpt;
  long nodata = poBand->GetNoDataValue();
  long pos;
  int iposx, iposy, iposy1;
  int State;
  float ix, ix2;
  // if (progress != NULL) progress->BeginProgress(CComBSTR("拷贝数据"));
  for (int i = 0; i < Height; i++) {
    lX = LeftTop.getX() - semiCellSize;
    for (int j = 0; j < Width; j++) {
      lX += CellSize;
      pX[j] = lX;
      pY[j] = lY;
    }
    if (poCT != NULL) poCT->Transform(Width, pX, pY);
    lX = LeftTop.getX() - semiCellSize;
    for (int j = 0; j < Width; j++) {
      lX += CellSize;

      if (!((PaintExtent.MinX < lX && PaintExtent.MaxX > lX)
            && (PaintExtent.MinY < lY && PaintExtent.MaxY > lY))) {
        Pos++;
        continue;
      }
      dpt.setX((pX[j] - FileExtent.MinX) / CellXSize);
      dpt.setY((pY[j] - FileExtent.MaxY) / CellYSize);
      if ((dpt.getX() < 0) || (dpt.getX() >= pW) || (dpt.getY() < 0) || (dpt.getY() >= pH)) {
        Pos++;
        continue;
      }
      iposx = dpt.getX();
      iposy = dpt.getY();
      if (dpt.getY() - iposy < 0.5) {
        iposy1 = iposy - 1;
        if (iposy1 < 0) iposy1 = iposy;
      } else {
        iposy1 = iposy + 1;
        if (iposy1 >= pH) iposy1 = iposy;
      }
      pos = iposy * pW + iposx;
      float v0;
      data->GetValueAsFloat(pos, &v0);
      if ((long)v0 == nodata)
        ix = nodata;
      else {
        if (dpt.getX() - iposx < 0.5) {
          if (iposx > 0) {
            float v;
            data->GetValueAsFloat(pos - 1, &v);
            if ((long)v != nodata)
              State = 1;
            else
              State = 2;
          } else
            State = 2;
          if (State == 2) {
            if (iposx + 1 >= pW)
              State = 0;
            else {
              float v;
              data->GetValueAsFloat(pos + 1, &v);
              if ((long)v == nodata) State = 0;
            }
          }
        } else {
          if (iposx + 1 < pW) {
            float v;
            data->GetValueAsFloat(pos + 1, &v);
            if ((long)v != nodata)
              State = 2;
            else
              State = 1;
          } else
            State = 1;
          if (State == 1) {
            if (iposx < 1)
              State = 0;
            else {
              float v;
              data->GetValueAsFloat(pos - 1, &v);
              if ((long)v == nodata) State = 0;
            }
          }
        }
        switch (State) {
          case 0: {
            ix = v0;
            break;
          }
          case 1: {
            float v;
            data->GetValueAsFloat(pos - 1, &v);
            ix = (dpt.getX() - iposx + 0.5) * (v0 - v) + v;
            break;
          }
          case 2: {
            float v;
            data->GetValueAsFloat(pos + 1, &v);
            ix = (dpt.getX() - iposx - 0.5) * (v - v0) + v0;
            break;
          }
        }
      }
      if (iposy1 == iposy) {
        pData->SetValueAsFloat(Pos, ix);
        Pos++;
        continue;
      }
      pos = iposy1 * pW + iposx;
      data->GetValueAsFloat(pos, &v0);
      if ((long)v0 == nodata)
        ix2 = nodata;
      else {
        if (dpt.getX() - iposx < 0.5) {
          if (iposx > 0) {
            float v;
            data->GetValueAsFloat(pos - 1, &v);
            if ((long)v != nodata)
              State = 1;
            else
              State = 2;
          } else
            State = 2;
          if (State == 2) {
            if (iposx + 1 >= pW)
              State = 0;
            else {
              float v;
              data->GetValueAsFloat(pos + 1, &v);
              if ((long)v == nodata) State = 0;
            }
          }
        } else {
          if (iposx + 1 < pW) {
            float v;
            data->GetValueAsFloat(pos + 1, &v);
            if ((long)v != nodata)
              State = 2;
            else
              State = 1;
          } else
            State = 1;
          if (State == 1) {
            if (iposx < 1)
              State = 0;
            else {
              float v;
              data->GetValueAsFloat(pos - 1, &v);
              if ((long)v == nodata) State = 0;
            }
          }
        }
        switch (State) {
          case 0: {
            ix2 = v0;
            break;
          }
          case 1: {
            float v;
            data->GetValueAsFloat(pos - 1, &v);
            ix2 = (dpt.getX() - iposx + 0.5) * (v0 - v) + v;
            break;
          }
          case 2: {
            float v;
            data->GetValueAsFloat(pos + 1, &v);
            ix2 = (dpt.getX() - iposx - 0.5) * (v - v0) + v0;
            break;
          }
        }
      }
      State = 0;
      if ((long)ix == nodata) State = 1;
      if ((long)ix2 == nodata) State += 2;
      switch (State) {
        case 0:
          pData->SetValueAsFloat(Pos,
                                 (ix2 - ix) * (dpt.getY() - iposy - 0.5) / (iposy1 - iposy) + ix);
          break;
        case 1:
          pData->SetValueAsFloat(Pos, ix2);
          break;
        case 2:
          pData->SetValueAsFloat(Pos, ix);
          break;
        case 3:
          pData->SetValueAsFloat(Pos, NoData);
          break;
      }
      Pos++;
    }
    lY -= CellSize;
    // if (progress != NULL) progress->SetPos((float)i / Height * 100);
  }
  // data->Release();
  delete[] pX;
  delete[] pY;
  return true;
}

/*--------------------------IFileFloatArray数据类型的定义----------------------------------*/
IFileFloatArray::~IFileFloatArray() {
  if (!this->rasterDataArray) {
    delete[] this->rasterDataArray;
  }
}
// IFileFloatArray::IFileFloatArray(void) {
//   this->bufferSize = fabsBuffer32MB;
//   this->size = 0;
// }
void IFileFloatArray::GetSize(long* pVal) { *pVal = this->size; }
void IFileFloatArray::SetSize(long size, float_t initialValue, bool* pVal) {
  this->size = size;
  // 动态分配数组大小
  this->rasterDataArray = new float_t[this->size]{initialValue};
  *pVal = true;
}

void IFileFloatArray::SetDefaultValue(float_t initialValue) {
  for (long i = 0; i < this->size; i++) {
    this->rasterDataArray[i] = initialValue;
  }
}

void IFileFloatArray::GetValueAsFloat(long pos, float_t* pVal) {
  *pVal = this->rasterDataArray[pos];
}

void IFileFloatArray::SetValueAsFloat(long pos, float_t newVal) {
  this->rasterDataArray[pos] = newVal;
}

void IFileFloatArray::GetBufferSize(FileArrayBufferSize* pVal) { *pVal = this->bufferSize; }
void IFileFloatArray::SetBufferSize(FileArrayBufferSize newVal) { this->bufferSize = newVal; }

float_t* IFileFloatArray::GetRasterDataArray() { return this->rasterDataArray; }
