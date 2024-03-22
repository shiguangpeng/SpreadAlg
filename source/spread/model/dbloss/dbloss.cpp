/**
 * @copyright ALL COPYRIGH RESERVED BY SHIGP
 * @author shigp
 */

#include <spatdata/spatdata.h>
#include <spread/model/combine.h>

using dbloss::CDBLossElement;

float CDBLossElement::GetZValue(double X, double Y) {
  int col = (X - XMin) / CellSize;
  if (col < 0) {
    return 0;
  } else if (col >= Cols) {
    return 0;
  }
  int row = (YMax - Y) / CellSize;
  if (row < 0)
    return 0;
  else if (row >= Rows)
    return 0;
  float z;
  pElevData->GetValueAsFloat(row * Cols + col, &z);
  if ((int64_t)z == (int64_t)NoData) {
    z = 0;
  }
  return z;
}

// void CDBLossElement::PrepareAnalyseEnvi() {
//   pAnalyse->get_AnalyseEnvi(&pEnvi);
//   pEnvi->get_Cols(&Cols);
//   pEnvi->get_Rows(&Rows);
//   IPoint* ppt;
//   pEnvi->get_LeftTop(&ppt);
//   ppt->get_X(&XMin);
//   ppt->get_Y(&YMax);
//   pEnvi->get_CellSize(&CellSize);
//   ppt->Release();
//   pAnalyse->get_NoData(&NoData);
//   pAnalyse->get_ElevationData(&pElevData);
//   pAnalyse->get_hm(&hm);
// }
