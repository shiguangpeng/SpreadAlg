/**
 * @copyright ALL COPYRIGH RESERVED BY SHIGP
 * @author shigp
 */
#include <spatdata/spatdata.h>
#include <spread/model/spreadanalyse.h>

using spatdata::CGDALRasterReaderByFileArray;
using spatdata::IGDALRasterProperties;
using spatdata::PBand_T;
using spread::spreadanalyse::AnalyseEnvironment;
using spread::spreadanalyse::CSpreadAnalyse;

AnalyseEnvironment::AnalyseEnvironment() {
  // 初始化类的成员变量
  leftTop = new OGRPoint;
  cellSize = 0.0;
  cols = rows = 0;
  outPutCellSize = 0.0;
}
// 重写接口对成员属性的get/set方法，对成员变量获取或设置值
void AnalyseEnvironment::GetLeftTop(OGRPoint **point) { *point = leftTop; }
void AnalyseEnvironment::SetLeftTop(OGRPoint *point) {
  leftTop->setX(point->getX());
  leftTop->setY(point->getY());
  leftTop->setZ(point->getZ());
}

void AnalyseEnvironment::GetCellSize(double_t *cellSize) {
  *cellSize = this->cellSize;
}
void AnalyseEnvironment::SetCellSize(double_t cellSize) {
  this->cellSize = cellSize;
}

void AnalyseEnvironment::GetCols(int64_t *cols) { *cols = this->cols; }
void AnalyseEnvironment::SetCols(int64_t cols) { this->cols = cols; }

void AnalyseEnvironment::GetRows(int64_t *rows) { *rows = this->rows; }
void AnalyseEnvironment::SetRows(int64_t rows) { this->rows = rows; }

void AnalyseEnvironment::GetOutputCellSize(double_t *outPutCellSize) {
  *outPutCellSize = this->outPutCellSize;
}
void AnalyseEnvironment::SetOutputCellSize(double_t outPutCellSize) {
  this->outPutCellSize = outPutCellSize;
}

CSpreadAnalyse::CSpreadAnalyse() {
  pEnvi = nullptr;
  errorInfo = "";
  // 初始化pElevs对象
  pElevs = nullptr;
  pElevData = nullptr;
  isInit = false;
  otherDataPreload = false;
}

/**
 * @brief 所有场强分析算法的基类
 * @return
 */
bool CSpreadAnalyse::InitEnvironment() {
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

  IGDALRasterProperties *pPro =
      dynamic_cast<CGDALRasterReaderByFileArray *>(this->pElevs);
  pPro->GetNoData(&noData);

  pPro->GetCols(&cols);
  pPro->GetRows(&rows);
  OGREnvelope *Extent;
  pPro->GetExtent(&Extent);
  OGRPoint *ppt = new OGRPoint;
  double_t left, top;
  left = Extent->MinX;
  top = Extent->MaxY;
  ppt->setX(left);
  ppt->setY(top);
  pPro->GetCellSize(&cellSize);
  if (pEnvi == nullptr) {
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
    OGRPoint *lt;
    pEnvi->GetLeftTop(&lt);
    if (lt != nullptr) {
      ppt = lt;
    } else {
      pEnvi->SetLeftTop(ppt);
    }
    int64_t rRows;
    int64_t rCols;
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
  // TODO(shigp): 在堆中分配空间
  pElevData = new IFileFloatArray;
  OGRPoint *lt;
  pEnvi->GetLeftTop(&lt);
  pElevs->GetBlockDataByCoord(lt, cellSize, cols, rows, noData, &pElevData);
  isInit = true;
  return true;
}
