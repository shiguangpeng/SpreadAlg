/**
 * @copyright ALL COPYRIGH RESERVED BY SHIGP
 * @author shigp
 */

#include <spatdata/spatdata.h>

#include <iostream>
#include <string>

using datatype::Station;
using spatdata::CGDALRasterReaderByFileArray;
using spatdata::CGDALRasterReaderByPixel;
using spatdata::CRect;
using spatdata::CStations;
using spatdata::IFileFloatArray;
using std::cout;
using std::endl;

void CStations::AddStation(Station* station) {
  this->stations.push_back(station);
}
void CStations::GetCount(int64_t* pVal) { *pVal = this->stations.size(); }
void CStations::RemoveAt(int64_t index) {
  // 获取数组开始的迭代器
  std::vector<Station*>::iterator spec = this->stations.begin() + index;
  this->stations.erase(spec);
}

void CStations::RemoveAll(void) { this->stations.clear(); }

void CStations::GetItem(int64_t index, Station* pVal) {
  // std::vector<spread::Station*>::iterator spec = this->stations.begin() +
  // index;
  *pVal = *this->stations.at(index);
}

void CStations::SetItem(int64_t index, Station newVal) {
  std::vector<Station*>::iterator spec = this->stations.begin() + index;
  this->stations.insert(spec, &newVal);
}

void CStations::PickHeightFromDEM(IGDALRasterReaderByPixel* pReader,
                                  bool* pVal) {
  int Size = stations.size();
  IGDALRasterProperties* pPro = new CGDALRasterReaderByFileArray;
  double CellSize;
  pPro->GetCellSize(&CellSize);
  OGREnvelope* rect;
  pPro->GetExtent(&rect);

  double XMin, YMax;
  XMin = rect->MinX;
  YMax = rect->MaxY;
  int64_t rows, cols;
  pPro->GetRows(&rows);
  pPro->GetCols(&cols);

  // 获取指定栅格的dem值
  float_t dV;
  for (int k = 0; k < Size; k++) {
    Station* para = stations.at(k);
    int Col = (para->x - XMin) / CellSize;
    int Row = (YMax - para->y) / CellSize;
    if ((Col < 0) || (Col >= cols) || (Row < 0) || (Row >= rows)) {
      para->dem = 0;
    } else {
      pReader->GetPixelValue(Col, Row, &dV);
      para->dem = dV;
    }
  }
  *pVal = true;
}

CRect::CRect() = default;
CRect::CRect(int left, int top, int right, int bottom) {
  this->left = left;
  this->top = top;
  this->bottom = bottom;
  this->right = right;
}
int CRect::Width() { return this->right - this->left; }
int CRect::Height() { return this->bottom - this->top; }

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
  }
  return pColorTable;
}
GDALRasterBand* CGDALRasterReaderByPixel::GetGDALBand() { return this->poBand; }
void CGDALRasterReaderByPixel::OpenRaster(std::string lpszPathName,
                                          bool* pVal) {
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
  poDataset = reinterpret_cast<GDALDataset*>(
      GDALOpen(lpszPathName.c_str(), GA_ReadOnly));
  if (poDataset == nullptr) {
    *pVal = false;
    return;
  }

  // GDALGetDriverShortName((GDALDriverH)poDataset); The SUBDATASETS domain
  // holds a list of child datasets. Normally this is used to provide pointers
  // to a list of images stored within a single multi image file.
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
void CGDALRasterReaderByPixel::GetPixelValue(int64_t col, int64_t row,
                                             float_t* data) {
  if (poBand == nullptr) {
    *data = -INT_MAX;
    return;
  }
  CPLErr errorNumber =
      poBand->RasterIO(GF_Read, col, row, 1, 1, data, 1, 1, GDT_Float32, 0, 0);
  cout << "获取结果标志：" << errorNumber << std::endl;
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

    poSubDataset =
        reinterpret_cast<GDALDataset*>(GDALOpen(pBand.rasterPath, GA_ReadOnly));
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
void CGDALRasterReaderByPixel::GetPathName(std::string* pVal) {
  *pVal = this->lpszPathName;
}
void CGDALRasterReaderByPixel::GetCurrentBand(PBand_T* pVal) { *pVal = pBand; }
void CGDALRasterReaderByPixel::GetRows(int64_t* pVal) { *pVal = this->rows; }
void CGDALRasterReaderByPixel::GetCols(int64_t* pVal) { *pVal = this->cols; }
void CGDALRasterReaderByPixel::GetExtent(OGREnvelope** pVal) {
  *pVal = this->extent;
}
void CGDALRasterReaderByPixel::GetCellSize(double_t* pVal) {
  *pVal = this->cellSize;
}
void CGDALRasterReaderByPixel::GetBandCount(int64_t* pVal) {
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
void CGDALRasterReaderByPixel::GetMinMax(bool bApproxOK, double_t* min,
                                         double_t* max, bool* pVal) {
  if (poBand == nullptr) {
    *pVal = false;
    return;
  }

  double_t mean, stddev;
  CPLErr pErr =
      poBand->GetStatistics(bApproxOK, true, min, max, &mean, &stddev);
  if (pErr == CE_None)
    *pVal = true;
  else
    *pVal = false;
  return;
}
void CGDALRasterReaderByPixel::ComputeRasterMinMax(bool bApproxOK,
                                                   double_t* min, double_t* max,
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
  } else {
    *pVal = false;
  }
  return;
}
void CGDALRasterReaderByPixel::ComputeStatistics(bool bApproxOK, double_t* min,
                                                 double_t* max, double_t* mean,
                                                 double_t* stddev, bool* pVal) {
  if (poBand == nullptr) {
    *pVal = false;
    return;
  }
  CPLErr pErr = poBand->ComputeStatistics(bApproxOK, min, max, mean, stddev,
                                          nullptr, nullptr);
  if (pErr == CE_None)
    *pVal = true;
  else
    *pVal = false;
  return;
}

void CGDALRasterReaderByPixel::GetStatistics(bool bApproxOK, double_t* min,
                                             double_t* max, double_t* mean,
                                             double_t* stddev, bool* pVal) {
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
  } else {
    pDataset = poDataset;
  }
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
  for (int64_t k = 0; k < static_cast<int>(metadata.size()); k++) {
    std::string sV = metadata.at(k);
    (*pVal)->push_back(sV);
  }
  return;
}

// TODO(shigp): 未实现
void CGDALRasterReaderByPixel::GetColorTable(
    std::vector<GDALColorTable>** pVal) {
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
  for (int64_t k = 0; k < size; k++) {
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

GDALRasterBand* CGDALRasterReaderByFileArray::GetGDALBand() {
  return this->poBand;
}

// todo: 空实现
OGREnvelope CGDALRasterReaderByFileArray::TransformRect(
    OGRCoordinateTransformation* poCT, OGREnvelope rt) {
  if (poCT == NULL) {
    return rt;
  }
  return rt;
}
OGREnvelope CGDALRasterReaderByFileArray::MapToPixelCoord(
    OGREnvelope MapExtent) {
  OGREnvelope PixelRect;
  PixelRect.MinX = (MapExtent.MinX - XMin) / CellSize;
  PixelRect.MinY = (MapExtent.MaxX - XMin) / CellSize;
  PixelRect.MaxY = (YMax - MapExtent.MaxY) / CellSize;
  PixelRect.MinY = (YMax - MapExtent.MinY) / CellSize;
  return PixelRect;
}
OGREnvelope CGDALRasterReaderByFileArray::PixelToMapCoord(
    OGREnvelope PixelExtent) {
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

void CGDALRasterReaderByFileArray::OpenRaster(std::string lpszPathName,
                                              bool* pVal) {
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
  poDataset = reinterpret_cast<GDALDataset*>(
      GDALOpen(lpszPathName.c_str(), GA_ReadOnly));
  if (poDataset == NULL) {
    *pVal = false;
    return;
  }

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
    if (poSubDataset != nullptr) {
      delete poSubDataset;
    }
    poSubDataset = nullptr;
    poSubDataset =
        reinterpret_cast<GDALDataset*>(GDALOpen(pBand.rasterPath, GA_ReadOnly));
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

void CGDALRasterReaderByFileArray::GetPathName(std::string* pVal) {
  *pVal = lpszPathName;
}

void CGDALRasterReaderByFileArray::GetCurrentBand(PBand_T* pVal) {
  *pVal = pBand;
}

void CGDALRasterReaderByFileArray::GetRows(int64_t* pVal) { *pVal = rows; }
void CGDALRasterReaderByFileArray::GetCols(int64_t* pVal) { *pVal = cols; }

// TODO(shigp): 返回了局部变量的引用，需要清理
void CGDALRasterReaderByFileArray::GetExtent(OGREnvelope** pVal) {
  OGREnvelope* pNew = new OGREnvelope();
  pNew->MaxX = XMax;
  pNew->MaxY = YMax;
  pNew->MinX = XMin;
  pNew->MinY = YMin;
  *pVal = pNew;
  return;
}
void CGDALRasterReaderByFileArray::GetCellSize(double_t* pVal) {
  *pVal = this->CellSize;
}
void CGDALRasterReaderByFileArray::GetBandCount(int64_t* pVal) {
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

void CGDALRasterReaderByFileArray::GetMinMax(bool bApproxOK, double_t* min,
                                             double_t* max, bool* pVal) {
  if (poBand == NULL) {
    *pVal = false;
    return;
  }
  double_t mean, stddev;
  CPLErr pErr =
      poBand->GetStatistics(bApproxOK, true, min, max, &mean, &stddev);
  if (pErr == CE_None)
    *pVal = true;
  else
    *pVal = false;
  return;
}

void CGDALRasterReaderByFileArray::GetStatistics(bool bApproxOK, double_t* min,
                                                 double_t* max, double_t* mean,
                                                 double_t* stddev, bool* pVal) {
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

void CGDALRasterReaderByFileArray::ComputeRasterMinMax(bool bApproxOK,
                                                       double_t* min,
                                                       double_t* max,
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
  } else {
    *pVal = false;
  }
  return;
}

void CGDALRasterReaderByFileArray::ComputeStatistics(
    bool bApproxOK, double_t* min, double_t* max, double_t* mean,
    double_t* stddev, bool* pVal) {
  if (poBand == NULL) {
    *pVal = false;
    return;
  }
  CPLErr pErr =
      poBand->ComputeStatistics(bApproxOK, min, max, mean, stddev, NULL, NULL);
  if (pErr == CE_None) {
    *pVal = true;
  } else {
    *pVal = false;
  }
  return;
}

// TODO(shigp): 清理堆区资源
void CGDALRasterReaderByFileArray::GetSpatialReference(
    OGRSpatialReference** pVal) {
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
    // TODO(shigp): 有返回值
    IsOk = (*pVal)->importFromWkt(info);
    if (IsOk != 0) {
      return;
    }
    return;
  }
}
void CGDALRasterReaderByFileArray::GetMetaData(
    std::vector<std::string>** pVal) {
  if (metadatas.size() == 0L) {
    *pVal = NULL;
    return;
  }
  std::vector<std::string> psa;
  // int64_t index;
  for (int64_t k = 0; k < (int64_t)metadatas.size(); k++) {
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

// TODO(shigp): 未实现GetColorTable
void CGDALRasterReaderByFileArray::GetColorTable(
    std::vector<GDALColorTable>** pVal) {
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
  // for (int64_t k = 0; k < Size; k++) {
  //   GDALColorEntry gce;
  //   gct->GetColorEntryAsRGB(k, &gce);

  //   (*pVal)->push_back(gct);
  // }
  return;
}

void CGDALRasterReaderByFileArray::GetBlockData(int64_t x1, int64_t y1,
                                                int64_t x2, int64_t y2,
                                                int64_t buffx, int64_t buffy,
                                                IFileFloatArray** pVal) {
  if (poBand == NULL) {
    *pVal = NULL;
    return;
  }
  if (x1 > x2) {
    int64_t temp = x1;
    x1 = x2;
    x2 = temp;
  }
  if (y1 > y2) {
    int64_t temp = y1;
    y1 = y2;
    y2 = temp;
  }
  int64_t buffer = buffx * buffy;
  if (buffer <= 0) {
    *pVal = nullptr;
    return;
  }
  if (buffer != BufferSize) {
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
  float_t difx = (x2 - x1 + 1) / buffx;
  float_t dify = (y2 - y1 + 1) / buffy;
  float_t* fV = new float_t[x2 - x1 + 1];
  int64_t CurrentRow = y1 - 1;
  int64_t pos = 0;
  for (float_t row = y1 + dify / 2; row < y2 + 1; row += dify) {
    if (CurrentRow != (int64_t)(row)) {
      CPLErr err = poBand->RasterIO(GF_Read, x1, row, (x2 - x1 + 1), 1, fV,
                                    (x2 - x1 + 1), 1, GDT_Float32, 0, 0);
      if (err != CE_None) {
        return;
      }
    }
    for (float_t col = x1 + difx / 2; col < x2 + 1; col += difx) {
      pData->SetValueAsFloat(pos, fV[static_cast<int64_t>(col)]);
      pos++;
    }
    CurrentRow = row;
  }
  *pVal = pData;
  delete[] fV;
  return;
}

void CGDALRasterReaderByFileArray::GetInterpolatedBlockData(
    int64_t x1, int64_t y1, int64_t x2, int64_t y2, int64_t BuffX,
    int64_t BuffY, IFileFloatArray** pVal) {
  int Width = BuffX;
  int Height = BuffY;
  CRect ImageRect;
  int64_t nodata = poBand->GetNoDataValue();
  int64_t buffer = Width * Height;
  if (buffer != BufferSize) {
    if (pData != NULL) {
      pData = NULL;
    }

    bool IsOk;
    pData->SetSize(buffer, nodata, &IsOk);
    if (!IsOk) {
      pData = NULL;
      BufferSize = 0;
      return;
    }
    BufferSize = buffer;
  }
  ImageRect = CRect(x1, y1, x2, y2);
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
    if (!ReadDataBlock(ImageRect.left, ImageRect.top, ImageRect.right,
                       ImageRect.bottom, Width, Height, pData)) {
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
    IFileFloatArray* data = GetDataBlock(
        ImageRect.left, ImageRect.top, ImageRect.right, ImageRect.bottom,
        ImageRect.Width() + 1, ImageRect.Height() + 1);
    if (data == NULL) {
      pData->SetDefaultValue(nodata);
      return;
    }
    int64_t Pos = 0;
    float_t posx, posy;
    int iposx, iposy, iposy1;
    int State;  // 0--None;1--Left;2--Right;
    int64_t Posi;
    int W = ImageRect.Width();
    int H = ImageRect.Height();
    float_t ix, ix2;
    float_t ratiox = (W + 1) / Width;
    float_t ratioy = (H + 1) / Height;
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
        float_t v0;
        data->GetValueAsFloat(Posi, &v0);
        if ((int64_t)v0 == nodata) {
          ix = nodata;
        } else {
          if (posx - iposx < 0.5) {
            if (iposx > 0) {
              float_t v;
              data->GetValueAsFloat(Posi - 1, &v);
              if ((int64_t)v != nodata) {
                State = 1;
              } else {
                State = 2;
              }
            } else {
              State = 2;
            }
            if (State == 2) {
              if (iposx + 1 > W) {
                State = 0;
              } else {
                float_t v;
                data->GetValueAsFloat(Posi + 1, &v);
                if ((int64_t)v == nodata) {
                  State = 0;
                }
              }
            }
          } else {
            if (iposx + 1 <= W) {
              float_t v;
              data->GetValueAsFloat(Posi + 1, &v);
              if ((int64_t)v != nodata) {
                State = 2;
              } else {
                State = 1;
              }
            } else {
              State = 1;
            }
            if (State == 1) {
              if (iposx < 1) {
                State = 0;
              } else {
                float_t v;
                data->GetValueAsFloat(Posi - 1, &v);
                if ((int64_t)v == nodata) {
                  State = 0;
                }
              }
            }
          }
          switch (State) {
            case 0: {
              ix = v0;
              break;
            }
            case 1: {
              float_t v;
              data->GetValueAsFloat(Posi - 1, &v);
              ix = (posx - iposx + 0.5) * (v0 - v) + v;
              break;
            }
            case 2: {
              float_t v;
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
        if ((int64_t)v0 == nodata) {
          ix2 = nodata;
        } else {
          if (posx - iposx < 0.5) {
            if (iposx > 0) {
              float_t v;
              data->GetValueAsFloat(Posi - 1, &v);
              if ((int64_t)v != nodata) {
                State = 1;
              } else {
                State = 2;
              }
            } else {
              State = 2;
            }
            if (State == 2) {
              if (iposx + 1 > W) {
                State = 0;
              } else {
                float_t v;
                data->GetValueAsFloat(Posi + 1, &v);
                if ((int64_t)v == nodata) State = 0;
              }
            }
          } else {
            if (iposx + 1 <= W) {
              float_t v;
              data->GetValueAsFloat(Posi + 1, &v);
              if ((int64_t)v != nodata)
                State = 2;
              else
                State = 1;
            } else {
              State = 1;
            }
            if (State == 1) {
              if (iposx < 1) {
                State = 0;
              } else {
                float_t v;
                data->GetValueAsFloat(Posi - 1, &v);
                if ((int64_t)v == nodata) State = 0;
              }
            }
          }
          switch (State) {
            case 0: {
              ix2 = v0;
              break;
            }
            case 1: {
              float_t v;
              data->GetValueAsFloat(Posi - 1, &v);
              ix2 = (posx - iposx + 0.5) * (v0 - v) + v;
              break;
            }
            case 2: {
              float_t v;
              data->GetValueAsFloat(Posi + 1, &v);
              ix2 = (posx - iposx - 0.5) * (v - v0) + v0;
              break;
            }
          }
        }
        State = 0;
        if ((int64_t)ix == nodata) State = 1;
        if ((int64_t)ix2 == nodata) State += 2;
        switch (State) {
          case 0:
            pData->SetValueAsFloat(
                Pos, (ix2 - ix) * (posy - iposy - 0.5) / (iposy1 - iposy) + ix);
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
    }
  }
  *pVal = pData;
  return;
}

void CGDALRasterReaderByFileArray::SetRefTargetSpatialReference(
    OGRSpatialReference* newVal) {
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
    {
      sp2.importFromWkt(wkt);
    }
  } else {
    return;
  }
  poCT = OGRCreateCoordinateTransformation(&sp2, &sp1);
  tpoCT = OGRCreateCoordinateTransformation(&sp1, &sp2);
  return;
}

void CGDALRasterReaderByFileArray::GetBlockDataByCoord(
    OGRPoint* LeftTop, double_t CellSize, int64_t Width, int64_t Height,
    float_t NoData, IFileFloatArray** pVal) {
  OGRPoint dpt;
  dpt.setX(LeftTop->getX());
  dpt.setY(LeftTop->getY());
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

void CGDALRasterReaderByFileArray::GetInterpolatedBlockDataByCoord(
    OGRPoint* LeftTop, double_t CellSize, int64_t Width, int64_t Height,
    float_t NoData, IFileFloatArray** pVal) {
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

bool CGDALRasterReaderByFileArray::GetDataBlock(OGRPoint LeftTop,
                                                float_t cellSize, int Width,
                                                int Height, float_t NoData) {
  OGREnvelope FullExtent;
  FullExtent.MinX = XMin;
  FullExtent.MinY = YMin;
  FullExtent.MaxX = XMax;
  FullExtent.MaxY = YMax;
  // DRect(XMin, YMax, XMax, YMin);
  float_t semiCellSize = cellSize / 2;
  OGREnvelope CurrentExtent;
  CurrentExtent.MinX = LeftTop.getX() + semiCellSize;
  CurrentExtent.MinY = LeftTop.getY() - Height * cellSize + semiCellSize;
  CurrentExtent.MaxX = LeftTop.getX() + cellSize * Width - semiCellSize;
  CurrentExtent.MaxY = LeftTop.getY() - semiCellSize;
  // (LeftTop.X + semiCellSize, LeftTop.Y - semiCellSize, LeftTop.X + cellSize
  // * Width - semiCellSize,
  //  LeftTop.Y - Height * cellSize + semiCellSize);
  CRect ImageRect;
  CRect TargetRect;
  int64_t nodata = poBand->GetNoDataValue();
  int64_t buffer = Width * Height;
  if (buffer != BufferSize) {
    if (pData != nullptr) {
      pData = nullptr;
    }
    BufferSize = buffer;
  }
  // FIXME(shigp): 在堆上动态分配空间记得在析构函数中清理
  pData = new IFileFloatArray;
  bool flag = false;
  pData->SetSize(BufferSize, 0, &flag);
  float_t* pvData = this->pData->GetRasterDataArray();
  if (!flag) {
    return false;
  }

  if (FullExtent.Contains(CurrentExtent)) {
    ImageRect.left = (CurrentExtent.MinX - XMin) / CellSize;
    ImageRect.top = (YMax - CurrentExtent.MaxY) / CellSize;
    ImageRect.right = (CurrentExtent.MaxX - XMin) / CellSize;
    ImageRect.bottom = (YMax - CurrentExtent.MinY) / CellSize;
    if ((ImageRect.Width() + 1 >= Width) ||
        (ImageRect.Height() + 1 >= Height)) {
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top,
                           ImageRect.Width() + 1, ImageRect.Height() + 1,
                           pvData, Width, Height, GDT_Float32, 0, 0) != CE_None)
        return false;
      if (nodata != (int64_t)NoData) {
        int64_t Pos = 0;
        for (int i = 0; i < Height; i++) {
          for (int j = 0; j < Width; j++) {
            if ((int64_t)pvData[Pos] == nodata) {
              pvData[Pos] = NoData;
            }
            Pos++;
          }
        }
      }
    } else {
      float_t* data = reinterpret_cast<float_t*>(
          CPLMalloc(sizeof(float_t) * (ImageRect.Width() + 1) *
                    sizeof(float_t) * (ImageRect.Height() + 1)));
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top,
                           ImageRect.Width() + 1, ImageRect.Height() + 1, data,
                           ImageRect.Width() + 1, ImageRect.Height() + 1,
                           GDT_Float32, 0, 0) != CE_None) {
        CPLFree(data);
        int64_t Size = Width * Height;
        for (int64_t k = 0; k < Size; k++) pvData[k] = NoData;
        return false;
      }
      OGREnvelope togr;
      togr.MinX = ImageRect.left;
      togr.MaxY = ImageRect.top;
      togr.MaxX = ImageRect.right;
      togr.MinY = ImageRect.bottom;
      OGREnvelope ext = PixelToMapCoord(togr);
      float_t Y = LeftTop.getY() - semiCellSize;
      int64_t Pos = 0;
      float_t X;
      int posx, posy;
      int64_t Posi;
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
          if ((int64_t)data[Posi] == nodata)
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
    int64_t Size = Width * Height;
    for (int64_t k = 0; k < Size; k++) {
      pvData[k] = NoData;
    }
    if (!FullExtent.Intersects(CurrentExtent)) return true;
    // TODO(shigp): 谁与谁相交
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
    if ((ImageRect.Width() >= TargetRect.Width()) ||
        (ImageRect.Height() >= TargetRect.Height())) {
      float_t* data = reinterpret_cast<float_t*>(
          CPLMalloc(sizeof(float_t) * (TargetRect.Width() + 1) *
                    sizeof(float_t) * (TargetRect.Height() + 1)));
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top,
                           ImageRect.Width() + 1, ImageRect.Height() + 1, data,
                           TargetRect.Width() + 1, TargetRect.Height() + 1,
                           GDT_Float32, 0, 0) != CE_None) {
        CPLFree(data);
        return false;
      }
      int64_t Pos;
      int64_t Posi = 0;
      for (int i = TargetRect.top; i <= TargetRect.bottom; i++) {
        Pos = i * Width + TargetRect.left;
        for (int j = TargetRect.left; j <= TargetRect.right; j++) {
          if ((int64_t)data[Posi] == nodata)
            pvData[Pos] = NoData;
          else
            pvData[Pos] = data[Posi];
          Pos++;
          Posi++;
        }
      }
      CPLFree(data);
    } else {
      float_t* data = reinterpret_cast<float_t*>(
          CPLMalloc(sizeof(float_t) * (ImageRect.Width() + 1) *
                    sizeof(float_t) * (ImageRect.Height() + 1)));
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top,
                           ImageRect.Width() + 1, ImageRect.Height() + 1, data,
                           ImageRect.Width() + 1, ImageRect.Height() + 1,
                           GDT_Float32, 0, 0) != CE_None) {
        CPLFree(data);
        return false;
      }
      int64_t Pos;
      OGREnvelope tmp;
      tmp.MinX = ImageRect.left;
      tmp.MaxY = ImageRect.top;
      tmp.MaxX = ImageRect.right;
      tmp.MinY = ImageRect.bottom;
      OGREnvelope ext = PixelToMapCoord(tmp);
      float_t Y = LeftTop.getY() - semiCellSize - TargetRect.top * cellSize;
      float_t X;
      int posx, posy;
      int64_t Posi;
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
          if ((int64_t)data[Posi] == nodata)
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

bool CGDALRasterReaderByFileArray::ReadDataBlock(int64_t x1, int64_t y1,
                                                 int64_t x2, int64_t y2,
                                                 int64_t buffx, int64_t buffy,
                                                 IFileFloatArray* pData) {
  if (poBand == nullptr) {
    return false;
  }
  if (x1 > x2) {
    int64_t temp = x1;
    x1 = x2;
    x2 = temp;
  }
  if (y1 > y2) {
    int64_t temp = y1;
    y1 = y2;
    y2 = temp;
  }
  int64_t buffer = buffx * buffy;
  if (buffer <= 0) {
    return false;
  }
  float_t difx = (float_t)(x2 - x1 + 1) / buffx;
  float_t dify = (float_t)(y2 - y1 + 1) / buffy;
  float_t* fV = new float_t[x2 - x1 + 1];
  int64_t CurrentRow = y1 - 1;
  int64_t pos = 0;
  for (float_t row = y1 + dify / 2; row < y2 + 1; row += dify) {
    if (CurrentRow != static_cast<int64_t>(row)) {
      CPLErr error = poBand->RasterIO(GF_Read, x1, row, (x2 - x1 + 1), 1, fV,
                                      (x2 - x1 + 1), 1, GDT_Float32, 0, 0);
      if (error != CE_None) {
        return false;
      }
    }
    for (float_t col = x1 + difx / 2; col < x2 + 1; col += difx) {
      pData->SetValueAsFloat(pos, fV[static_cast<int64_t>(col) - x1]);
      pos++;
    }
    CurrentRow = row;
  }
  delete[] fV;
  return true;
}

IFileFloatArray* CGDALRasterReaderByFileArray::GetDataBlock(
    int64_t x1, int64_t y1, int64_t x2, int64_t y2, int64_t buffx,
    int64_t buffy) {
  if (poBand == NULL) {
    return NULL;
  }
  if (x1 > x2) {
    int64_t temp = x1;
    x1 = x2;
    x2 = temp;
  }
  if (y1 > y2) {
    int64_t temp = y1;
    y1 = y2;
    y2 = temp;
  }
  int64_t buffer = buffx * buffy;
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
  float_t difx = (float_t)(x2 - x1 + 1) / buffx;
  float_t dify = (float_t)(y2 - y1 + 1) / buffy;
  float_t* fV = new float_t[x2 - x1 + 1];
  int64_t CurrentRow = y1 - 1;
  int64_t pos = 0;
  for (float_t row = y1 + dify / 2; row < y2 + 1; row += dify) {
    if (CurrentRow != static_cast<int64_t>(row)) {
      CPLErr err = poBand->RasterIO(GF_Read, x1, row, (x2 - x1 + 1), 1, fV,
                                    (x2 - x1 + 1), 1, GDT_Float32, 0, 0);
      if (err != CE_None) {
        return nullptr;
      }
    }
    for (float_t col = x1 + difx / 2; col < x2 + 1; col += difx) {
      pData->SetValueAsFloat(pos, fV[static_cast<int64_t>(col) - x1]);
      pos++;
    }
    CurrentRow = row;
  }
  delete[] fV;
  return pData;
}

bool CGDALRasterReaderByFileArray::GetInterpolatedDataBlock(
    OGRPoint LeftTop, float_t cellSize, int Width, int Height, float_t NoData) {
  // DRect FullExtent(XMin, YMax, XMax, YMin);
  OGREnvelope FullExtent;
  FullExtent.MinX = XMin;
  FullExtent.MinY = YMin;
  FullExtent.MaxX = XMax;
  FullExtent.MaxY = YMax;
  float_t semiCellSize = cellSize / 2;

  OGREnvelope CurrentExtent;
  CurrentExtent.MinX = LeftTop.getX() + semiCellSize;
  CurrentExtent.MaxY = LeftTop.getY() - semiCellSize;
  CurrentExtent.MaxX = LeftTop.getX() + cellSize * Width - semiCellSize;
  CurrentExtent.MinY = LeftTop.getY() - Height * cellSize + semiCellSize;

  CRect ImageRect;
  CRect TargetRect;
  int64_t nodata = poBand->GetNoDataValue();
  int64_t buffer = Width * Height;
  if (buffer != BufferSize) {
    if (pData != NULL) {
      pData = nullptr;
    }
    BufferSize = buffer;
  }
  // 判断完成后，需要创建IFileFloatArray类型的实例
  // FIXME(shigp): 在堆上动态分配空间记得在析构函数中清理
  pData = new IFileFloatArray;
  bool flag = false;
  // float_t* pvData = (float_t*)pData->pvData;
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
    if ((ImageRect.Width() + 1 >= Width) ||
        (ImageRect.Height() + 1 >= Height)) {
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top,
                           ImageRect.Width() + 1, ImageRect.Height() + 1,
                           pvData, Width, Height, GDT_Float32, 0,
                           0) != CE_None) {
        return false;
      }
      if ((int64_t)nodata != (int64_t)NoData) {
        int64_t Pos = 0;
        for (int i = 0; i < Height; i++) {
          for (int j = 0; j < Width; j++) {
            if ((int64_t)pvData[Pos] == nodata) pvData[Pos] = NoData;
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
      float_t* data = reinterpret_cast<float_t*>(
          CPLMalloc(sizeof(float_t) * (ImageRect.Width() + 1) *
                    sizeof(float_t) * (ImageRect.Height() + 1)));
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top,
                           ImageRect.Width() + 1, ImageRect.Height() + 1, data,
                           ImageRect.Width() + 1, ImageRect.Height() + 1,
                           GDT_Float32, 0, 0) != CE_None) {
        CPLFree(data);
        int64_t Size = Width * Height;
        for (int64_t k = 0; k < Size; k++) pvData[k] = NoData;
        return false;
      }

      OGREnvelope envelope;
      envelope.MinX = ImageRect.left;
      envelope.MinX = ImageRect.left;
      envelope.MinX = ImageRect.left;
      envelope.MinX = ImageRect.left;
      OGREnvelope ext = PixelToMapCoord(envelope);
      float_t Y = LeftTop.getY() - semiCellSize;
      int64_t Pos = 0;
      float_t X;
      float_t posx, posy;
      int iposx, iposy, iposy1;
      int State;  // 0--None;1--Left;2--Right;
      int64_t Posi;
      int W = ImageRect.Width();
      int H = ImageRect.Height();
      float_t ix, ix2;
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
          if ((int64_t)data[Posi] == nodata) {
            ix = nodata;
          } else {
            if (posx - iposx < 0.5) {
              if (iposx > 0) {
                if ((int64_t)data[Posi - 1] != nodata) {
                  State = 1;
                } else {
                  State = 2;
                }
              } else {
                State = 2;
              }
              if (State == 2) {
                if (iposx + 1 > W)
                  State = 0;
                else if ((int64_t)data[Posi + 1] == nodata)
                  State = 0;
              }
            } else {
              if (iposx + 1 <= W) {
                if ((int64_t)data[Posi + 1] != nodata) {
                  State = 2;
                } else {
                  State = 1;
                }
              } else {
                State = 1;
              }
              if (State == 1) {
                if (iposx < 1)
                  State = 0;
                else if ((int64_t)data[Posi - 1] == nodata)
                  State = 0;
              }
            }
            switch (State) {
              case 0: {
                ix = data[Posi];
                break;
              }
              case 1: {
                ix = (posx - iposx + 0.5) * (data[Posi] - data[Posi - 1]) +
                     data[Posi - 1];
                break;
              }
              case 2: {
                ix = (posx - iposx - 0.5) * (data[Posi + 1] - data[Posi]) +
                     data[Posi];
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
          if ((int64_t)data[Posi] == nodata) {
            ix2 = nodata;
          } else {
            if (posx - iposx < 0.5) {
              if (iposx > 0) {
                if ((int64_t)data[Posi - 1] != nodata) {
                  State = 1;
                } else {
                  State = 2;
                }
              } else {
                State = 2;
              }
              if (State == 2) {
                if (iposx + 1 > W)
                  State = 0;
                else if ((int64_t)data[Posi + 1] == nodata)
                  State = 0;
              }
            } else {
              if (iposx + 1 <= W) {
                if ((int64_t)data[Posi + 1] != nodata) {
                  State = 2;
                } else {
                  State = 1;
                }
              } else {
                State = 1;
              }
              if (State == 1) {
                if (iposx < 1)
                  State = 0;
                else if ((int64_t)data[Posi - 1] == nodata)
                  State = 0;
              }
            }
            switch (State) {
              case 0: {
                ix2 = data[Posi];
                break;
              }
              case 1: {
                ix2 = (posx - iposx + 0.5) * (data[Posi] - data[Posi - 1]) +
                      data[Posi - 1];
                break;
              }
              case 2: {
                ix2 = (posx - iposx - 0.5) * (data[Posi + 1] - data[Posi]) +
                      data[Posi];
                break;
              }
            }
          }
          State = 0;
          if ((int64_t)ix == nodata) State = 1;
          if ((int64_t)ix2 == nodata) State += 2;
          switch (State) {
            case 0:
              pvData[Pos] =
                  (ix2 - ix) * (posy - iposy - 0.5) / (iposy1 - iposy) + ix;
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
    int64_t Size = Width * Height;
    for (int64_t k = 0; k < Size; k++) pvData[k] = NoData;
    // TODO(shigp): 判断前一个envelope是否与后面个有交集
    if (!FullExtent.Intersects(CurrentExtent)) return true;
    // TODO(shigp): 取交集
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
    if ((ImageRect.Width() >= TargetRect.Width()) ||
        (ImageRect.Height() >= TargetRect.Height())) {
      float_t* data = reinterpret_cast<float_t*>(
          CPLMalloc(sizeof(float_t) * (TargetRect.Width() + 1) *
                    sizeof(float_t) * (TargetRect.Height() + 1)));
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top,
                           ImageRect.Width() + 1, ImageRect.Height() + 1, data,
                           TargetRect.Width() + 1, TargetRect.Height() + 1,
                           GDT_Float32, 0, 0) != CE_None) {
        CPLFree(data);
        return false;
      }
      int64_t Pos;
      int64_t Posi = 0;
      for (int i = TargetRect.top; i <= TargetRect.bottom; i++) {
        Pos = i * Width + TargetRect.left;
        for (int j = TargetRect.left; j <= TargetRect.right; j++) {
          if ((int64_t)data[Posi] == (int64_t)nodata)
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
      // TODO(shigp): ??这是什么意思
      ImageRect.left -= 1;
      ImageRect.top -= 1;
      ImageRect.right += 1;
      ImageRect.bottom += 1;
      if (ImageRect.left < 0) ImageRect.left++;
      if (ImageRect.right >= cols) ImageRect.right--;
      if (ImageRect.top < 0) ImageRect.top++;
      if (ImageRect.bottom >= rows) ImageRect.bottom--;
      float_t* data = static_cast<float_t*>(
          CPLMalloc(sizeof(float_t) * (ImageRect.Width() + 1) *
                    sizeof(float_t) * (ImageRect.Height() + 1)));
      if (poBand->RasterIO(GF_Read, ImageRect.left, ImageRect.top,
                           ImageRect.Width() + 1, ImageRect.Height() + 1, data,
                           ImageRect.Width() + 1, ImageRect.Height() + 1,
                           GDT_Float32, 0, 0) != CE_None) {
        CPLFree(data);
        return false;
      }
      int64_t Pos;
      OGREnvelope ogrEnvelope;
      ogrEnvelope.MinX = ImageRect.left;
      ogrEnvelope.MinY = ImageRect.bottom;
      ogrEnvelope.MaxX = ImageRect.right;
      ogrEnvelope.MaxY = ImageRect.top;
      // OGREnvelope ext = PixelToMapCoord(ImageRect);
      OGREnvelope ext = PixelToMapCoord(ogrEnvelope);
      float_t Y = LeftTop.getY() - semiCellSize - TargetRect.top * cellSize;
      float_t X;
      float_t posx, posy;
      int iposx, iposy, iposy1;
      int State;  // 0--None;1--Left;2--Right;
      int64_t Posi;
      int W = cols;
      int H = rows;
      float_t ix, ix2;
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
          if ((int64_t)data[Posi] == nodata) {
            ix = nodata;
          } else {
            if (posx - iposx < 0.5) {
              if (iposx > 0) {
                if ((int64_t)data[Posi - 1] != nodata) {
                  State = 1;
                } else {
                  State = 2;
                }
              } else {
                State = 2;
              }
              if (State == 2) {
                if (iposx + 1 > W)
                  State = 0;
                else if ((int64_t)data[Posi + 1] == nodata)
                  State = 0;
              }
            } else {
              if (iposx + 1 <= W) {
                if ((int64_t)data[Posi + 1] != nodata)
                  State = 2;
                else
                  State = 1;
              } else {
                State = 1;
              }
              if (State == 1) {
                if (iposx < 1)
                  State = 0;
                else if ((int64_t)data[Posi - 1] == nodata)
                  State = 0;
              }
            }
            switch (State) {
              case 0: {
                ix = data[Posi];
                break;
              }
              case 1: {
                ix = (posx - iposx + 0.5) * (data[Posi] - data[Posi - 1]) +
                     data[Posi - 1];
                break;
              }
              case 2: {
                ix = (posx - iposx - 0.5) * (data[Posi + 1] - data[Posi]) +
                     data[Posi];
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
          if ((int64_t)data[Posi] == nodata) {
            ix2 = nodata;
          } else {
            if (posx - iposx < 0.5) {
              if (iposx > 0) {
                if ((int64_t)data[Posi - 1] != nodata) {
                  State = 1;
                } else {
                  State = 2;
                }
              } else {
                State = 2;
              }
              if (State == 2) {
                if (iposx + 1 > W)
                  State = 0;
                else if ((int64_t)data[Posi + 1] == nodata)
                  State = 0;
              }
            } else {
              if (iposx + 1 <= W) {
                if ((int64_t)data[Posi + 1] != nodata) {
                  State = 2;
                } else {
                  State = 1;
                }
              } else {
                State = 1;
              }
              if (State == 1) {
                if (iposx < 1)
                  State = 0;
                else if ((int64_t)data[Posi - 1] == nodata)
                  State = 0;
              }
            }
            switch (State) {
              case 0: {
                ix2 = data[Posi];
                break;
              }
              case 1: {
                ix2 = (posx - iposx + 0.5) * (data[Posi] - data[Posi - 1]) +
                      data[Posi - 1];
                break;
              }
              case 2: {
                ix2 = (posx - iposx - 0.5) * (data[Posi + 1] - data[Posi]) +
                      data[Posi];
                break;
              }
            }
          }
          State = 0;
          if ((int64_t)ix == nodata) State = 1;
          if ((int64_t)ix2 == nodata) State += 2;
          switch (State) {
            case 0:
              pvData[Pos] =
                  (ix2 - ix) * (posy - iposy - 0.5) / (iposy1 - iposy) + ix;
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

bool CGDALRasterReaderByFileArray::GetDataBlockWithProj(OGRPoint LeftTop,
                                                        float_t CellSize,
                                                        int Width, int Height,
                                                        float_t NoData) {
  OGREnvelope MapExtent;
  MapExtent.MinX = LeftTop.getX();
  MapExtent.MaxY = LeftTop.getY();
  MapExtent.MaxX = LeftTop.getX() + CellSize * Width;
  MapExtent.MinY = LeftTop.getY() - Height * CellSize;

  OGREnvelope FileExtent =
      TransformRect(poCT, MapExtent);  // 将地图范围转换为图层范围
  OGREnvelope pixelRect =
      MapToPixelCoord(FileExtent);  // 由图层范围得到图层像素范围
  CRect PixelRect;
  PixelRect.left = pixelRect.MinX;
  PixelRect.top = pixelRect.MaxY;
  PixelRect.right = pixelRect.MaxX;
  PixelRect.bottom = pixelRect.MinY;
  // PixelRect.NormalizeRect();
  // TODO(shigp): 归一化矩形？
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
  // int64_t Size = Width * Height;
  int64_t buffer = Width * Height;
  if (buffer != BufferSize) {
    if (pData != NULL) {
      // pData->Release();
      pData = NULL;
    }
    bool IsOk;
    pData->SetSize(buffer, NoData, &IsOk);
    if (!IsOk) {
      pData = NULL;
      BufferSize = 0;
      return false;
    }
    BufferSize = buffer;
  } else {
    pData->SetDefaultValue(NoData);
  }
  // 图层像素范围进行纠正
  OGREnvelope envelope;
  envelope.MinX = PixelRect.left;
  envelope.MaxY = PixelRect.top;
  envelope.MaxX = PixelRect.right;
  envelope.MinY = PixelRect.bottom;
  FileExtent = PixelToMapCoord(envelope);
  // 由图层像素范围反算图层范围
  OGREnvelope PaintExtent =
      TransformRect(tpoCT, FileExtent);  // 将图层范围转换为要显示的地图范围
  int pW, pH;  // 表示文件中与PixelRect成比例的宽和高，依据W和H计算
  float_t r1 = (float_t)Width / Height;
  float_t r2 = (float_t)(PixelRect.Width() + 1) / (PixelRect.Height() + 1);
  if (r1 >= r2) {
    if (Width < PixelRect.Width() + 1) {
      pW = Width;
      float_t fpH = pW / r2;
      pH = fpH;
      if (fpH - pH >= 0.5) pH++;
    } else {
      pW = PixelRect.Width() + 1;
      pH = PixelRect.Height() + 1;
    }
  } else {
    if (Height < PixelRect.Height() + 1) {
      pH = Height;
      float_t fpW = pH * r2;
      pW = fpW;
      if (fpW - pW >= 0.5) pW++;
    } else {
      pW = PixelRect.Width() + 1;
      pH = PixelRect.Height() + 1;
    }
  }
  IFileFloatArray* data = GetDataBlock(
      PixelRect.left, PixelRect.top, PixelRect.right, PixelRect.bottom, pW, pH);
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
  float_t semiCellSize = CellSize / 2;
  int64_t Pos = 0;
  double CellXSize = (FileExtent.MaxX - FileExtent.MinX) / pW;
  double CellYSize = (FileExtent.MaxY - FileExtent.MinY) / pH;
  lY = LeftTop.getY() - semiCellSize;
  OGRPoint ppt;
  int64_t nodata = poBand->GetNoDataValue();
  int64_t pos;
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
      if (!(((PaintExtent.MinX < lX) && (PaintExtent.MaxX > lX)) &&
            ((PaintExtent.MinY < lY) && (PaintExtent.MaxY > lY)))) {
        Pos++;
        continue;
      }
      ppt.setX((pX[j] - FileExtent.MinX) / CellXSize);
      ppt.setY((pY[j] - FileExtent.MaxY) / CellYSize);
      if ((ppt.getX() < 0) || (ppt.getX() >= pW) || (ppt.getY() < 0) ||
          (ppt.getY() >= pH)) {
        Pos++;
        continue;
      } else {
        pos = ppt.getY() * pW + ppt.getX();
        float_t v;
        data->GetValueAsFloat(pos, &v);
        if ((int64_t)v == nodata)
          pData->SetValueAsFloat(Pos, NoData);
        else
          pData->SetValueAsFloat(Pos, v);
      }
      Pos++;
    }
    lY -= CellSize;
    // if (progress != NULL) progress->SetPos((float_t)i / Height * 100);
  }
  // data->Release();
  delete[] pX;
  delete[] pY;
  return true;
}
bool CGDALRasterReaderByFileArray::GetInterpolatedDataBlockWithProj(
    OGRPoint LeftTop, float_t CellSize, int Width, int Height, float_t NoData) {
  OGREnvelope MapExtent;
  MapExtent.MinX = LeftTop.getX();
  MapExtent.MaxX = LeftTop.getX() + CellSize * Width;
  MapExtent.MinY = LeftTop.getY() - Height * CellSize;
  MapExtent.MaxY = LeftTop.getY();

  OGREnvelope FileExtent =
      TransformRect(poCT, MapExtent);  // 将地图范围转换为图层范围
  OGREnvelope pixelRect =
      MapToPixelCoord(FileExtent);  // 由图层范围得到图层像素范围
  CRect PixelRect;
  PixelRect.left = pixelRect.MinX;
  PixelRect.top = pixelRect.MaxY;
  PixelRect.right = pixelRect.MaxX;
  PixelRect.bottom = pixelRect.MinY;
  // PixelRect.NormalizeRect();
  // PixelRect.InflateRect(1, 1, 1, 1);
  // TODO(shigp): 归一化矩形？
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
  // int64_t Size = Width * Height;
  int64_t buffer = Width * Height;
  if (PixelRect.left > PixelRect.right) return false;
  if (PixelRect.top > PixelRect.bottom) return false;
  if (buffer != BufferSize) {
    if (pData != NULL) {
      // pData->Release();
      pData = NULL;
    }
    // ::CoCreateInstance(CLSID_FileFloatArray, NULL, CLSCTX_INPROC_SERVER,
    // IID_IFileFloatArray,
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
  } else {
    pData->SetDefaultValue(NoData);
  }
  // 图层像素范围进行纠正
  OGREnvelope envelope;
  envelope.MinX = PixelRect.left;
  envelope.MinY = PixelRect.bottom;
  envelope.MaxX = PixelRect.right;
  envelope.MaxY = PixelRect.top;
  FileExtent = PixelToMapCoord(envelope);
  // 由图层像素范围反算图层范围
  OGREnvelope PaintExtent =
      TransformRect(tpoCT, FileExtent);  // 将图层范围转换为要显示的地图范围
  int pW, pH;  // 表示文件中与PixelRect成比例的宽和高，依据W和H计算
  float_t r1 = (float_t)Width / Height;
  float_t r2 = (float_t)(PixelRect.Width() + 1) / (PixelRect.Height() + 1);
  if (r1 >= r2) {
    if (Width < PixelRect.Width() + 1) {
      pW = Width;
      float_t fpH = pW / r2;
      pH = fpH;
      if (fpH - pH >= 0.5) pH++;
    } else {
      pW = PixelRect.Width() + 1;
      pH = PixelRect.Height() + 1;
    }
  } else {
    if (Height < PixelRect.Height() + 1) {
      pH = Height;
      float_t fpW = pH * r2;
      pW = fpW;
      if (fpW - pW >= 0.5) pW++;
    } else {
      pW = PixelRect.Width() + 1;
      pH = PixelRect.Height() + 1;
    }
  }
  IFileFloatArray* data = GetDataBlock(
      PixelRect.left, PixelRect.top, PixelRect.right, PixelRect.bottom, pW, pH);
  if (data == nullptr) {
    pData = nullptr;
    BufferSize = 0;
    return false;
  }
  double *pX, *pY;
  pX = new double[Width];
  pY = new double[Width];
  double lX, lY;
  float_t semiCellSize = CellSize / 2;
  int64_t Pos = 0;
  double CellXSize = (FileExtent.MaxX - FileExtent.MinX) / pW;
  double CellYSize = (FileExtent.MaxY - FileExtent.MinY) / pH;
  lY = LeftTop.getY() - semiCellSize;
  OGRPoint dpt;
  int64_t nodata = poBand->GetNoDataValue();
  int64_t pos;
  int iposx, iposy, iposy1;
  int State;
  float_t ix, ix2;
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

      if (!((PaintExtent.MinX < lX && PaintExtent.MaxX > lX) &&
            (PaintExtent.MinY < lY && PaintExtent.MaxY > lY))) {
        Pos++;
        continue;
      }
      dpt.setX((pX[j] - FileExtent.MinX) / CellXSize);
      dpt.setY((pY[j] - FileExtent.MaxY) / CellYSize);
      if ((dpt.getX() < 0) || (dpt.getX() >= pW) || (dpt.getY() < 0) ||
          (dpt.getY() >= pH)) {
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
      float_t v0;
      data->GetValueAsFloat(pos, &v0);
      if ((int64_t)v0 == nodata) {
        ix = nodata;
      } else {
        if (dpt.getX() - iposx < 0.5) {
          if (iposx > 0) {
            float_t v;
            data->GetValueAsFloat(pos - 1, &v);
            if ((int64_t)v != nodata) {
              State = 1;
            } else {
              State = 2;
            }
          } else {
            State = 2;
          }
          if (State == 2) {
            if (iposx + 1 >= pW) {
              State = 0;
            } else {
              float_t v;
              data->GetValueAsFloat(pos + 1, &v);
              if ((int64_t)v == nodata) State = 0;
            }
          }
        } else {
          if (iposx + 1 < pW) {
            float_t v;
            data->GetValueAsFloat(pos + 1, &v);
            if ((int64_t)v != nodata) {
              State = 2;
            } else {
              State = 1;
            }
          } else {
            State = 1;
          }
          if (State == 1) {
            if (iposx < 1) {
              State = 0;
            } else {
              float_t v;
              data->GetValueAsFloat(pos - 1, &v);
              if ((int64_t)v == nodata) {
                State = 0;
              }
            }
          }
        }
        switch (State) {
          case 0: {
            ix = v0;
            break;
          }
          case 1: {
            float_t v;
            data->GetValueAsFloat(pos - 1, &v);
            ix = (dpt.getX() - iposx + 0.5) * (v0 - v) + v;
            break;
          }
          case 2: {
            float_t v;
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
      if ((int64_t)v0 == nodata) {
        ix2 = nodata;
      } else {
        if (dpt.getX() - iposx < 0.5) {
          if (iposx > 0) {
            float_t v;
            data->GetValueAsFloat(pos - 1, &v);
            if ((int64_t)v != nodata) {
              State = 1;
            } else {
              State = 2;
            }
          } else {
            State = 2;
          }
          if (State == 2) {
            if (iposx + 1 >= pW) {
              State = 0;
            } else {
              float_t v;
              data->GetValueAsFloat(pos + 1, &v);
              if ((int64_t)v == nodata) State = 0;
            }
          }
        } else {
          if (iposx + 1 < pW) {
            float_t v;
            data->GetValueAsFloat(pos + 1, &v);
            if ((int64_t)v != nodata)
              State = 2;
            else
              State = 1;
          } else {
            State = 1;
          }
          if (State == 1) {
            if (iposx < 1) {
              State = 0;
            } else {
              float_t v;
              data->GetValueAsFloat(pos - 1, &v);
              if ((int64_t)v == nodata) State = 0;
            }
          }
        }
        switch (State) {
          case 0: {
            ix2 = v0;
            break;
          }
          case 1: {
            float_t v;
            data->GetValueAsFloat(pos - 1, &v);
            ix2 = (dpt.getX() - iposx + 0.5) * (v0 - v) + v;
            break;
          }
          case 2: {
            float_t v;
            data->GetValueAsFloat(pos + 1, &v);
            ix2 = (dpt.getX() - iposx - 0.5) * (v - v0) + v0;
            break;
          }
        }
      }
      State = 0;
      if ((int64_t)ix == nodata) State = 1;
      if ((int64_t)ix2 == nodata) State += 2;
      switch (State) {
        case 0:
          pData->SetValueAsFloat(
              Pos,
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
  }
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

void IFileFloatArray::GetSize(int64_t* pVal) { *pVal = this->size; }
void IFileFloatArray::SetSize(int64_t size, float_t initialValue, bool* pVal) {
  this->size = size;
  // 动态分配数组大小
  this->rasterDataArray = new float_t[this->size]{initialValue};
  *pVal = true;
}

void IFileFloatArray::SetDefaultValue(float_t initialValue) {
  for (int64_t i = 0; i < this->size; i++) {
    this->rasterDataArray[i] = initialValue;
  }
}

void IFileFloatArray::GetValueAsFloat(int64_t pos, float_t* pVal) {
  *pVal = this->rasterDataArray[pos];
}

void IFileFloatArray::SetValueAsFloat(int64_t pos, float_t newVal) {
  this->rasterDataArray[pos] = newVal;
}

void IFileFloatArray::GetBufferSize(FileArrayBufferSize* pVal) {
  *pVal = this->bufferSize;
}
void IFileFloatArray::SetBufferSize(FileArrayBufferSize newVal) {
  this->bufferSize = newVal;
}

float_t* IFileFloatArray::GetRasterDataArray() { return this->rasterDataArray; }
