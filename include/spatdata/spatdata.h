/**
 * @copyright ALL COPYRIGH RESERVED BY SHIGP
 * @author shigp
 */

#ifndef INCLUDE_SPATDATA_SPATDATA_H_
#define INCLUDE_SPATDATA_SPATDATA_H_
#include <gdal_priv.h>
#include <spatdata/datatype.h>

#include <cmath>
#include <string>
#include <vector>

using datatype::FileArrayBufferSize;
using datatype::RasterDataType;
using datatype::Station;
namespace spatdata {
// 为CGDALRasterReaderByXXX系列的列中的variant类型做一个替换，替换成这个联合
union PBand_T {
  int bandNumber;
  char *rasterPath;
};

struct IDataArray {
 public:
  // virtual std::string GetName(std::string pVal) = 0;
  // virtual void SetName(std::string newVal) = 0;
  virtual ~IDataArray() = default;
  virtual void GetSize(int64_t *pVal) = 0;
  virtual void SetSize(int64_t size, float_t initialValue, bool *pVal) = 0;

  virtual void SetDefaultValue(float_t initialValue = 0) = 0;

  virtual void GetValueAsFloat(int64_t pos, float_t *pVal) = 0;
  virtual void SetValueAsFloat(int64_t pos, float_t newVal) = 0;
};

struct IFileDataArray : public IDataArray {
  virtual ~IFileDataArray() = default;

 public:
  virtual void GetBufferSize(FileArrayBufferSize *pVal) = 0;
  virtual void SetBufferSize(FileArrayBufferSize newVal) = 0;
};

struct IGDALRasterReader {
 public:
  virtual ~IGDALRasterReader() = default;

 public:
  virtual void OpenRaster(std::string lpszPathName, bool *pVal) = 0;
  virtual void SetRasterBand(PBand_T pBand, bool *pVal) = 0;
};

// TODO(shigp): ATTENTION:原项目这里是一个空的结构体
// 原项目的IFileFloatArray没有找到实现，这个类型是封装了的一维数组，
// 直接使用该声明，将其修改成实现类（结构体），实现其接口规定的方法
// 这样不必动项目中有关该类型的代码
struct IFileFloatArray : public IFileDataArray {
 public:
  ~IFileFloatArray();

 protected:
  // 栅格数据数组
  int64_t size;
  // 一维数组指针，动态分配
  float_t *rasterDataArray;
  FileArrayBufferSize bufferSize;

 public:
  // std::string GetName(std::string pVal) override;
  // void SetName(std::string newVal) override;

  void GetSize(int64_t *pVal) override;
  void SetSize(int64_t size, float_t initialValue, bool *pVal) override;

  void SetDefaultValue(float_t initialValue = 0) override;

  void GetValueAsFloat(int64_t pos, float_t *pVal) override;
  void SetValueAsFloat(int64_t pos, float_t newVal) override;

  void GetBufferSize(FileArrayBufferSize *pVal) override;
  void SetBufferSize(FileArrayBufferSize newVal) override;
  float_t *GetRasterDataArray();
};

struct IGDALRasterReaderByFileArray : public IGDALRasterReader {
 public:
  virtual ~IGDALRasterReaderByFileArray() = default;

 public:
  virtual void GetBlockData(int64_t x1, int64_t y1, int64_t x2, int64_t y2,
                            int64_t buffx, int64_t buffy,
                            IFileFloatArray **pVal) = 0;

  virtual void GetInterpolatedBlockData(int64_t x1, int64_t y1, int64_t x2,
                                        int64_t y2, int64_t BuffX,
                                        int64_t BuffY,
                                        IFileFloatArray **pVal) = 0;

  // 替换为OGRSpatialReference
  virtual void SetRefTargetSpatialReference(OGRSpatialReference *newVal) = 0;

  virtual void GetBlockDataByCoord(OGRPoint *leftTop, double_t cellSize,
                                   int64_t width, int64_t height,
                                   float_t noData, IFileFloatArray **pVal) = 0;

  virtual void GetInterpolatedBlockDataByCoord(OGRPoint *LeftTop,
                                               double_t CellSize, int64_t Width,
                                               int64_t Height, float_t NoData,
                                               IFileFloatArray **pVal) = 0;
};

struct IGDALRasterProperties : public IGDALRasterReader {
 public:
  virtual ~IGDALRasterProperties() = default;

 public:
  virtual void GetPathName(std::string *pVal) = 0;
  virtual void GetCurrentBand(PBand_T *pVal) = 0;
  virtual void GetRows(int64_t *pVal) = 0;
  virtual void GetCols(int64_t *pVal) = 0;
  virtual void GetExtent(OGREnvelope **pVal) = 0;
  virtual void GetCellSize(double_t *pVal) = 0;
  virtual void GetBandCount(int64_t *pVal) = 0;
  virtual void GetNoData(double_t *pVal) = 0;
  virtual void SetNoData(double_t pVal) = 0;
  virtual void GetMetaData(std::vector<std::string> **pVal) = 0;
  virtual void GetDataType(RasterDataType *pVal) = 0;
  virtual void GetMinMax(bool bApproxOK, double_t *min, double_t *max,
                         bool *pVal) = 0;
  virtual void GetStatistics(bool bApproxOK, double_t *min, double_t *max,
                             double_t *mean, double_t *stddev, bool *pVal) = 0;
  virtual void ComputeRasterMinMax(bool bApproxOK, double_t *min, double_t *max,
                                   bool *pVal) = 0;
  virtual void ComputeStatistics(bool bApproxOK, double_t *min, double_t *max,
                                 double_t *mean, double_t *stddev,
                                 bool *pVal) = 0;
  virtual void GetSpatialReference(OGRSpatialReference **pVal) = 0;
  virtual void GetColorTable(std::vector<GDALColorTable> **pVal) = 0;
};
// 按像素读取栅格，接口
struct IGDALRasterReaderByPixel : public IGDALRasterReader {
 public:
  virtual void GetPixelValue(int64_t col, int64_t row, float_t *pVal) = 0;
};
// 环境类的接口
struct IAnalyseEnvironment {
 protected:
  /// @brief 获取栅格左上角坐标
  /// @return void
  virtual void GetLeftTop(OGRPoint **p) = 0;
  /// @brief 设置栅格左上角坐标
  /// @param point OGR坐标对象
  virtual void SetLeftTop(OGRPoint *point) = 0;

  /// @brief 获取栅格分辨率
  /// @return 栅格分辨率
  virtual void GetCellSize(double_t *cellSize) = 0;
  /// @brief 设置栅格分辨率
  /// @param cellSize 栅格分辨率
  virtual void SetCellSize(double_t cellSize) = 0;

  /// @brief 获取栅格的列数
  /// @return 栅格列总数
  virtual void GetCols(int64_t *cols) = 0;
  /// @brief 设置栅格的列数
  /// @param cols 栅格列数
  virtual void SetCols(int64_t cols) = 0;

  /// @brief 获取栅格行数
  /// @return 栅格行数
  virtual void GetRows(int64_t *rows) = 0;
  /// @brief 设置栅格行数
  /// @param rows 栅格行数
  virtual void SetRows(int64_t rows) = 0;

  /// @brief 获取输出栅格的分辨率
  /// @return 输出栅格的分辨率
  virtual void GetOutputCellSize(double_t *cellSize) = 0;
  /// @brief 设置输出栅格的分辨率
  /// @param outPutCellSize 输出栅格分辨率
  virtual void SetOutputCellSize(double_t outPutCellSize) = 0;
};

struct IStations {
 public:
  virtual void AddStation(Station *station) = 0;
  virtual void GetCount(int64_t *pVal) = 0;
  virtual void RemoveAt(int64_t index) = 0;
  virtual void RemoveAll(void) = 0;
  virtual void GetItem(int64_t index, Station *pVal) = 0;
  virtual void SetItem(int64_t index, Station newVal) = 0;
  virtual void PickHeightFromDEM(IGDALRasterReaderByPixel *pReader,
                                 bool *pVal) = 0;
};

struct CStations : public IStations {
 protected:
  std::vector<Station *> stations;

 public:
  void AddStation(Station *station) override;
  void GetCount(int64_t *pVal) override;
  void RemoveAt(int64_t index) override;
  void RemoveAll(void) override;
  void GetItem(int64_t index, Station *pVal) override;
  void SetItem(int64_t index, Station newVal) override;
  void PickHeightFromDEM(IGDALRasterReaderByPixel *pReader,
                         bool *pVal) override;
};

// 自定义矩形结构体类型
struct CRect {
 public:
  int left;
  int top;
  int right;
  int bottom;

 public:
  CRect();
  CRect(int left, int top, int right, int bottom);

 public:
  int Width();
  int Height();
};

// IGDALRasterReaderByPixel继承了IGDALRasterProperties和IGDALRasterReaderByPixel
class CGDALRasterReaderByPixel : public IGDALRasterProperties,
                                 public IGDALRasterReaderByPixel {
 protected:
  int64_t rows;
  int64_t cols;
  OGREnvelope *extent;
  double_t cellSize;
  GDALDataset *poDataset;
  GDALDataset *poSubDataset;
  // 由 VARIANT类型 pBand
  union PBand_T pBand;

  GDALRasterBand *poBand;
  std::vector<std::string> metadata;
  std::string lpszPathName;

 public:
  GDALColorTable *CreateColorTable();
  GDALRasterBand *GetGDALBand();

 public:
  void OpenRaster(std::string lpszPathName, bool *pVal) override;
  void SetRasterBand(PBand_T pBand, bool *pVal) override;
  void GetPathName(std::string *pVal) override;
  void GetCurrentBand(PBand_T *pVal) override;
  void GetRows(int64_t *pVal) override;
  void GetCols(int64_t *pVal) override;
  void GetExtent(OGREnvelope **pVal) override;
  void GetCellSize(double_t *pVal) override;
  void GetBandCount(int64_t *pVal) override;
  void GetNoData(double_t *pVal) override;
  void SetNoData(double_t pVal) override;
  void GetMinMax(bool bApproxOK, double_t *min, double_t *max,
                 bool *pVal) override;
  void GetSpatialReference(OGRSpatialReference **pVal) override;
  void GetDataType(RasterDataType *pVal) override;
  void GetPixelValue(int64_t col, int64_t row, float_t *data) override;
  void ComputeRasterMinMax(bool bApproxOK, double_t *min, double_t *max,
                           bool *pVal) override;
  void GetStatistics(bool bApproxOK, double_t *min, double_t *max,
                     double_t *mean, double_t *stddev, bool *pVal) override;
  void ComputeStatistics(bool bApproxOK, double_t *min, double_t *max,
                         double_t *mean, double_t *stddev, bool *pVal) override;
  void GetMetaData(std::vector<std::string> **pVal) override;
  void GetColorTable(std::vector<GDALColorTable> **pVal) override;
};

class CGDALRasterReaderByFileArray : public IGDALRasterProperties,
                                     public IGDALRasterReaderByFileArray {
 public:
  ~CGDALRasterReaderByFileArray() = default;

 protected:
  int64_t rows;
  int64_t cols;
  double XMin, YMin, XMax, YMax;
  double_t CellSize;
  GDALDataset *poDataset;
  GDALDataset *poSubDataset;
  GDALRasterBand *poBand;
  PBand_T pBand;
  std::vector<std::string> metadatas;
  std::string lpszPathName;
  IFileFloatArray *pData;
  int64_t BufferSize;
  OGRCoordinateTransformation *poCT;
  OGRCoordinateTransformation *tpoCT;

 public:
  CGDALRasterReaderByFileArray();
  GDALColorTable *CreateColorTable();
  GDALRasterBand *GetGDALBand();

 protected:
  OGREnvelope TransformRect(OGRCoordinateTransformation *poCT, OGREnvelope rt);
  OGREnvelope MapToPixelCoord(OGREnvelope MapExtent);
  OGREnvelope PixelToMapCoord(OGREnvelope PixelExtent);

  IFileFloatArray *GetDataBlock(int64_t x1, int64_t y1, int64_t x2, int64_t y2,
                                int64_t buffx, int64_t buffy);
  bool GetDataBlock(OGRPoint LeftTop, float CellSize, int Width, int Height,
                    float NoData);
  bool ReadDataBlock(int64_t x1, int64_t y1, int64_t x2, int64_t y2,
                     int64_t buffx, int64_t buffy, IFileFloatArray *pData);
  bool GetInterpolatedDataBlock(OGRPoint LeftTop, float CellSize, int Width,
                                int Height, float NoData);
  bool GetDataBlockWithProj(OGRPoint LeftTop, float CellSize, int Width,
                            int Height, float NoData);
  bool GetInterpolatedDataBlockWithProj(OGRPoint LeftTop, float CellSize,
                                        int Width, int Height, float NoData);

 public:
  void OpenRaster(std::string lpszPathName, bool *pVal) override;
  void SetRasterBand(PBand_T band, bool *pVal) override;
  void GetPathName(std::string *pVal) override;
  void GetCurrentBand(PBand_T *pVal) override;
  void GetRows(int64_t *pVal) override;
  void GetCols(int64_t *pVal) override;
  void GetExtent(OGREnvelope **pVal) override;
  void GetCellSize(double_t *pVal) override;
  void GetBandCount(int64_t *pVal) override;
  void GetNoData(double_t *pVal) override;
  void SetNoData(double_t pVal) override;
  void GetMinMax(bool bApproxOK, double_t *min, double_t *max,
                 bool *pVal) override;
  void GetStatistics(bool bApproxOK, double_t *min, double_t *max,
                     double_t *mean, double_t *stddev, bool *pVal) override;
  void ComputeRasterMinMax(bool bApproxOK, double_t *min, double_t *max,
                           bool *pVal) override;
  void ComputeStatistics(bool bApproxOK, double_t *min, double_t *max,
                         double_t *mean, double_t *stddev, bool *pVal) override;
  void GetSpatialReference(OGRSpatialReference **pVal) override;
  void GetMetaData(std::vector<std::string> **pVal) override;
  void GetDataType(RasterDataType *pVal) override;
  void GetColorTable(std::vector<GDALColorTable> **pVal) override;

  void GetBlockData(int64_t x1, int64_t y1, int64_t x2, int64_t y2,
                    int64_t buffx, int64_t buffy,
                    IFileFloatArray **pVal) override;
  void GetInterpolatedBlockData(int64_t x1, int64_t y1, int64_t x2, int64_t y2,
                                int64_t BuffX, int64_t BuffY,
                                IFileFloatArray **pVal) override;
  void SetRefTargetSpatialReference(OGRSpatialReference *newVal) override;
  void GetBlockDataByCoord(OGRPoint *LeftTop, double_t CellSize, int64_t Width,
                           int64_t Height, float NoData,
                           IFileFloatArray **pVal) override;
  void GetInterpolatedBlockDataByCoord(OGRPoint *LeftTop, double_t CellSize,
                                       int64_t Width, int64_t Height,
                                       float_t NoData,
                                       IFileFloatArray **pVal) override;
};

}  // namespace spatdata

#endif  // INCLUDE_SPATDATA_SPATDATA_H_
