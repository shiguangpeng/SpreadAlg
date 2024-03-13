#ifndef SPREADALG_INCLUDE_SPREAD_
#define SPREADALG_INCLUDE_SPREAD_

#include <error.h>
#include <spread/version.h>

#include <string>
#include <vector>

#include "gdal_priv.h"

namespace spread {
  /// @brief 描述站的结构体
  typedef struct StationST {
    /// @brief 台站id
    long id;
    /// @brief 东西方向坐标
    double_t x;
    /// @brief 南北方向坐标
    double_t y;
    /// @brief 经度
    double_t lontitude;
    /// @brief 纬度
    double_t latitude;
    /// @brief 台站的灵敏度
    double_t freqThreshold;
    /// @brief 站的高度
    double_t stationHeight;
    /// @brief 频率
    double_t frequency;
    /// @brief 功率
    double_t power;
    /// @brief 站所处高程
    /// @attention
    /// 我国的高程基准是青岛港验潮站的长期观测资料推算出的黄海平均海面（72.260m作为零高程面）
    double_t dem;
    /// @brief 站的名字
    std::string name;
    /// @brief 站的唯一标识
    std::string guid;
    /// @brief 方位角，水平方向
    double_t azimuth;
    /// @brief 高度角，竖直方向
    double_t altitude;
    /// @brief 天线直径
    double_t antennaDiameter;
  } Station;
  /** --------------------------- CDEMAnalyse对应的数据的封装 ----------------------------------*/
  struct IDEMAnalyse {
  public:
    virtual void getAnalyseEnvi() = 0;
  };
  struct ICombineAnalyse : public IDEMAnalyse {};

  /* ------------------ DBLossElement结构体，父结构体，保存绕射传播的各种参数
     -------------------*/
  struct IDBLossElement {
  public:
    // 是否复杂模型
    virtual void getComplicatedModel(_Bool *isComplicated) = 0;
    // 设置CombineAnalyse
    virtual void setCombineAnalyse(ICombineAnalyse *combineAnalyse) = 0;
    virtual void prepareAnalyseEnvi(void) = 0;
    virtual double_t GetDBLoss(Station station, OGRPoint point) = 0;
    virtual double_t GetDBLossRev(Station station, OGRPoint point) = 0;
  };

  /* ------------------------------------- 枚举类型定义 ----------------------------------------*/
  typedef enum FileDataArrayTypeST {
    fdatUnknown = 0,
    fdatByte = 0x1,
    fdatShort = 0x2,
    fdatInteger = 0x3,
    fdatFloat = 0x4
  } FileDataArrayType;

  typedef enum FileArrayBufferSizeST {
    fabsBufferHalfMB = 0,
    fabsBuffer1MB = 0x1,
    fabsBuffer2MB = 0x2,
    fabsBuffer4MB = 0x3,
    fabsBuffer8MB = 0x4,
    fabsBuffer16MB = 0x5,
    fabsBuffer32MB = 0x6,
    fabsBuffer64MB = 0x7
  } FileArrayBufferSize;

  // 栅格数据类型，可选tiff，erdas，envi，pci
  typedef enum RasterCreateFileTypeST {
    rcftTiff = 0x0,
    rcftErdasImagine = 0x1,
    rcftPCIPix = 0x2,
    rcftEnviHdr = 0x3
  } RasterCreateFileType;

  typedef enum RasterCreateDataTypeST {
    rcdtUnknown = 0,
    rcdtByte = 0x1,
    rcdtUInt16 = 0x2,
    rcdtInt16 = 0x3,
    rcdtUInt32 = 0x4,
    rcdtInt32 = 0x5,
    rcdtFloat32 = 0x6,
    rcdtFloat64 = 0x7
  } RasterCreateDataType;

  typedef enum RasterDataTypeST {
    rdtUnknown = 0,
    rdtByte = 0x1,
    rdtUInt16 = 0x2,
    rdtInt16 = 0x3,
    rdtUInt32 = 0x4,
    rdtInt32 = 0x5,
    rdtFloat32 = 0x6,
    rdtFloat64 = 0x7,
    rdtCInt16 = 0x8,
    rdtCInt32 = 0x9,
    rdtCFloat32 = 0xa,
    rdtCFloat64 = 0xb,
    rdtTwoValue = 0xc,
    rdtFourValue = 0xd,
    rdtEightValue = 0xe,
    rdtSixteenValue = 0xf
  } RasterDataType;

  // 为CGDALRasterReaderByXXX系列的列中的variant类型做一个替换，替换成这个联合
  union PBand_T {
    int bandNumber;
    char *rasterPath;
  };
  /** --------------------------- 数据类型封装-------------------------------------------------*/
  struct IDataArray {
  public:
    // virtual std::string GetName(std::string pVal) = 0;
    // virtual void SetName(std::string newVal) = 0;
    virtual ~IDataArray() = default;
    virtual void GetSize(long *pVal) = 0;
    virtual void SetSize(long size, float_t initialValue, bool *pVal) = 0;

    virtual void SetDefaultValue(float_t initialValue = 0) = 0;

    virtual void GetValueAsFloat(long pos, float_t *pVal) = 0;
    virtual void SetValueAsFloat(long pos, float_t newVal) = 0;

    // virtual void GetValuesAsByte(long fromPos, std::vector<double_t> **refData) = 0;
    // virtual void SetValuesAsByte(long fromPos, std::vector<double_t> *refData) = 0;
    // virtual void GetValuesAsShort(long fromPos, std::vector<double_t> **refData) = 0;
    // virtual void SetValuesAsShort(long fromPos, std::vector<double_t> *refData) = 0;
    // virtual void GetValuesAsInteger(long fromPos, std::vector<double_t> **refData) = 0;
    // virtual void SetValuesAsInteger(long fromPos, std::vector<double_t> *refData) = 0;
    // virtual void GetValuesAsFloat(long fromPos, std::vector<double_t> **refData) = 0;
    // virtual void SetValuesAsFloat(long fromPos, std::vector<double_t> *refData) = 0;
  };

  struct IFileDataArray : public IDataArray {
    virtual ~IFileDataArray() = default;

  public:
    // virtual void GetType(FileDataArrayType *pVal) = 0;
    // virtual void GetDirectory(std::string pVal) = 0;
    // virtual void SetDirectory(std::string newVal) = 0;
    virtual void GetBufferSize(FileArrayBufferSize *pVal) = 0;
    virtual void SetBufferSize(FileArrayBufferSize newVal) = 0;
  };
  /* ----------------------------- GDAL读取栅格数据的封装 ---------------------------------------*/
  struct IGDALRasterReader {
  public:
    virtual ~IGDALRasterReader() = default;

  public:
    virtual void OpenRaster(std::string lpszPathName, bool *pVal) = 0;
    virtual void SetRasterBand(PBand_T pBand, bool *pVal) = 0;
  };
  // FIXME: ATTENTION:原项目这里是一个空的结构体
  // 原项目的IFileFloatArray没有找到实现，这个类型是封装了的一维数组，直接使用该声明，将其修改成实现类（结构体），实现其接口规定的方法
  // 这样不必动项目中有关该类型的代码
  struct IFileFloatArray : public IFileDataArray {
  public:
    ~IFileFloatArray();

  protected:
    // 栅格数据数组
    long size;
    // 一维数组指针，动态分配
    float_t *rasterDataArray;
    FileArrayBufferSize bufferSize;

  public:
    // std::string GetName(std::string pVal) override;
    // void SetName(std::string newVal) override;

    void GetSize(long *pVal) override;
    void SetSize(long size, float_t initialValue, bool *pVal) override;

    void SetDefaultValue(float_t initialValue = 0) override;

    void GetValueAsFloat(long pos, float_t *pVal) override;
    void SetValueAsFloat(long pos, float_t newVal) override;

    void GetBufferSize(FileArrayBufferSize *pVal) override;
    void SetBufferSize(FileArrayBufferSize newVal) override;
    float_t *GetRasterDataArray();
  };

  struct IGDALRasterReaderByFileArray : public IGDALRasterReader {
  public:
    virtual ~IGDALRasterReaderByFileArray() = default;

  public:
    virtual void GetBlockData(long x1, long y1, long x2, long y2, long buffx, long buffy,
                              IFileFloatArray **pVal)
        = 0;

    virtual void GetInterpolatedBlockData(long x1, long y1, long x2, long y2, long BuffX,
                                          long BuffY, IFileFloatArray **pVal)
        = 0;

    // 替换为OGRSpatialReference
    virtual void SetRefTargetSpatialReference(OGRSpatialReference *newVal) = 0;

    virtual void GetBlockDataByCoord(OGRPoint *leftTop, double_t cellSize, long width, long height,
                                     float_t noData, IFileFloatArray **pVal)
        = 0;

    virtual void GetInterpolatedBlockDataByCoord(OGRPoint *LeftTop, double_t CellSize, long Width,
                                                 long Height, float_t NoData,
                                                 IFileFloatArray **pVal)
        = 0;
  };

  struct IGDALRasterProperties : public IGDALRasterReader {
  public:
    virtual ~IGDALRasterProperties() = default;

  public:
    virtual void GetPathName(std::string *pVal) = 0;
    virtual void GetCurrentBand(PBand_T *pVal) = 0;
    virtual void GetRows(long *pVal) = 0;
    virtual void GetCols(long *pVal) = 0;
    virtual void GetExtent(OGREnvelope **pVal) = 0;
    virtual void GetCellSize(double_t *pVal) = 0;
    virtual void GetBandCount(long *pVal) = 0;
    virtual void GetNoData(double_t *pVal) = 0;
    virtual void SetNoData(double_t pVal) = 0;
    virtual void GetMetaData(std::vector<std::string> **pVal) = 0;
    virtual void GetDataType(RasterDataType *pVal) = 0;
    virtual void GetMinMax(bool bApproxOK, double_t *min, double_t *max, bool *pVal) = 0;
    virtual void GetStatistics(bool bApproxOK, double_t *min, double_t *max, double_t *mean,
                               double_t *stddev, bool *pVal)
        = 0;
    virtual void ComputeRasterMinMax(bool bApproxOK, double_t *min, double_t *max, bool *pVal) = 0;
    virtual void ComputeStatistics(bool bApproxOK, double_t *min, double_t *max, double_t *mean,
                                   double_t *stddev, bool *pVal)
        = 0;
    virtual void GetSpatialReference(OGRSpatialReference **pVal) = 0;
    virtual void GetColorTable(std::vector<GDALColorTable> **pVal) = 0;
  };
  // 按像素读取栅格，接口
  struct IGDALRasterReaderByPixel : public IGDALRasterReader {
  public:
    virtual void GetPixelValue(long col, long row, float_t *pVal) = 0;
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
    virtual void GetCols(long *cols) = 0;
    /// @brief 设置栅格的列数
    /// @param cols 栅格列数
    virtual void SetCols(long cols) = 0;

    /// @brief 获取栅格行数
    /// @return 栅格行数
    virtual void GetRows(long *rows) = 0;
    /// @brief 设置栅格行数
    /// @param rows 栅格行数
    virtual void SetRows(long rows) = 0;

    /// @brief 获取输出栅格的分辨率
    /// @return 输出栅格的分辨率
    virtual void GetOutputCellSize(double_t *cellSize) = 0;
    /// @brief 设置输出栅格的分辨率
    /// @param outPutCellSize 输出栅格分辨率
    virtual void SetOutputCellSize(double_t outPutCellSize) = 0;
  };
  /* ------------------------------Stations ------------------------------------ */
  struct IStations {
  public:
    virtual void AddStation(Station *station) = 0;
    virtual void GetCount(long *pVal) = 0;
    virtual void RemoveAt(long index) = 0;
    virtual void RemoveAll(void) = 0;
    virtual void GetItem(long index, Station *pVal) = 0;
    virtual void SetItem(long index, Station newVal) = 0;
    virtual void PickHeightFromDEM(IGDALRasterReaderByPixel *pReader, bool *pVal) = 0;
  };

  struct CStations : public IStations {
  protected:
    std::vector<Station *> stations;

  public:
    void AddStation(Station *station) override;
    void GetCount(long *pVal) override;
    void RemoveAt(long index) override;
    void RemoveAll(void) override;
    void GetItem(long index, Station *pVal) override;
    void SetItem(long index, Station newVal) override;
    void PickHeightFromDEM(IGDALRasterReaderByPixel *pReader, bool *pVal) override;
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

  /**
   * @brief 场强分析算法的基类
   */
  class CSpreadAnalyse {
  public:
    CSpreadAnalyse();
    virtual ~CSpreadAnalyse() = default;

  public:
    /**
     * @brief 初始化分析环境
     */
    bool InitEnvironment();

  public:
    /**
     * @brief 分析环境参数结构体，封装了对环境参数的get与set
     */
    struct AnalyseEnvironment *pEnvi;

    /**
     * @brief 高程栅格路径
     */
    std::string elevationPath;

    // 封装过的根据条件获取栅格数据的对象
    IGDALRasterReaderByFileArray *pElevs;

    /**
     * @brief 错误信息
     */
    std::string errorInfo;

    // 栅格数据对象
    IFileFloatArray *pElevData;

    /**
     * @brief 栅格列数
     */
    long cols;

    /**
     * @brief 栅格行数
     */
    long rows;

    /**
     * @brief 栅格大小
     */
    double_t cellSize;

    /**
     * @brief 最小经度位置（左上角）
     */
    double_t xMin;

    /**
     * @brief 最大纬度位置（左上角）
     */
    double_t yMax;

    /**
     * @brief 栅格无数据（或无有效值）的填充值
     */
    double_t noData;

    /**
     * @brief 环境是否已经被初始化
     * @note true（已经初始化）/false（未初始化）
     */
    bool isInit;

    /**
     * @brief 场强分析所需要的其他必要数据是否已经被加载
     * @note true（已经加载）/false（未加载）
     */
    bool otherDataPreload;
  };

  class CFieldStrengthAnalyse : public CSpreadAnalyse {
  public:
    IStations *pStations;
    float_t hm;
    std::string outputPath;
    IFileFloatArray *pData;
    float_t offsetDB;
    bool needComputeAll;
    std::vector<double_t> crsv;
    Station stationInfo;
    /// @brief 分析范围，默认是整个dem，每个子类有自己的get/putExtent方法，为该数组设置值
    double_t subExtent[4];

  public:
    CFieldStrengthAnalyse();
    virtual ~CFieldStrengthAnalyse() = default;
    // 场强分析算法的调用入口，每个场强分析类都有自己的场强分析调用入口，重定义。
    bool FieldStrengthAnalyse(std::string savePath, RasterCreateFileType type);

  protected:
    /* 子类必须实现的方法 */
    /// @brief 发射站覆盖分析
    /// @param stationInfo 台站信息
    /// @param rsv 动态数组，用于保存每个子类保存的数据，每次使用前需要留意数组中已经存在的数据
    /// @param point 待分析的点坐标
    /// @return 返回点对点之间场强分析的结果
    virtual float_t GetRadiuValue(Station &stationInfo, std::vector<double_t> &rsv, OGRPoint &point)
        = 0;

    /// @brief
    /// @param stationInfo 台站信息
    /// @param rsv
    /// @param point 待分析的点坐标
    /// @return
    virtual float_t GetRadiuValueRev(Station &stationInfo, std::vector<double_t> &rsv,
                                     OGRPoint &point)
        = 0;

    /// @brief
    /// @param stationInfo
    /// @param rsv 是一个数组，该数组里面装了台站的初始的发射功率以及其他相关信息
    virtual void PrepareReservedValues(Station &stationInfo, std::vector<double_t> &rsv) = 0;

    /// @brief 分析所需要的其他栅格数据
    /// @return 返回准备的状态
    virtual bool PrepareOtherData() = 0;

    /// @brief 获取模型的名称
    /// @return 返回模型的名称
    virtual std::string GetModelName() = 0;

  protected:
    void ComputeOneStation(Station &para);
    bool ComputeOneStationByExtent(Station &para, std::vector<double_t> &rsv);
    long int *GetSubRowCol(double xmin, double ymin, double xmax, double ymax);
    // bool CFieldStrengthAnalyse::SaveOutput(IFileFloatArray *pData, RasterCreateFileType type,
    //                                        RasterCreateDataType dtype);
  };

  /**
   * @brief
   * 组合分析类，是需要多种数据源以及绕射模型计算的类的父类，子类中明确了特定的绕射模型，该类中根据子类传入的rsv数组中的绕射
   * 模型调用对应的实现计算处绕射衰减
   * @attention 该类的子类模型均有绕射衰减的计算，但是绕射衰减的计算是可选的，默认值为0
   */
  class CCombineAnalyse : public CFieldStrengthAnalyse {
  public:
    CCombineAnalyse() = default;
    ~CCombineAnalyse() = default;

  public:
    bool FieldStrengthAnalyse(std::string savePath, RasterCreateFileType type);

  protected:
    float_t GetRadiuValue(Station &stationInfo, std::vector<double_t> &rsv,
                          OGRPoint &point) override;

    float_t GetRadiuValueRev(Station &para, std::vector<double_t> &rsv, OGRPoint &point) override;
    // void PrepareReservedValues(Station &stationInfo, std::vector<double_t> &rsv) override;
    bool PrepareOtherData() override;
    std::string GetModelName() override;

  protected:
    std::vector<IGDALRasterReaderByFileArray *> otherReaders;
    std::vector<IFileFloatArray *> otherDatas;
    std::vector<IDBLossElement *> els;
    std::vector<std::string> otherNames;
    std::vector<std::string> otherPaths;
  };

  /* -----------------------------自由空间传播模型，暂只实现必要方法------------------------------*/

  class CFreeSpaceAnalyse : public CCombineAnalyse {
  protected:
    float_t GetRadiuValue(Station &stationInfo, std::vector<double_t> &rsv, OGRPoint &point);
    void PrepareReservedValues(Station &para, std::vector<double_t> &rsv);

  public:
    // 调用父类CFreeSpaceAnalyse的FieldStrengthAnalyse方法
    void FieldStrengthAnalyse(std::string savePath, RasterCreateFileType type, bool *pVal);
  };

  /**
   * @brief 分析环境的参数结构体，专门管理分析环境参数
   */
  struct AnalyseEnvironment : public IAnalyseEnvironment {
  public:
    AnalyseEnvironment();
    // ~AnalyseEnvironment();

  private:
    OGRPoint *leftTop;
    double_t cellSize;
    long cols;
    long rows;
    double_t outPutCellSize;

  public:
    /// @brief 获取栅格左上角坐标
    /// @return OGRPoint*
    void GetLeftTop(OGRPoint **p) override;
    /// @brief 设置栅格左上角坐标
    /// @param point OGR坐标对象
    void SetLeftTop(OGRPoint *point) override;
    /// @brief 获取栅格分辨率
    /// @return 栅格分辨率
    void GetCellSize(double_t *cellSize) override;
    /// @brief 设置栅格分辨率
    /// @param cellSize 栅格分辨率
    void SetCellSize(double_t cellSize) override;
    /// @brief 获取栅格的列数
    /// @return 栅格列总数
    void GetCols(long *cols) override;
    /// @brief 设置栅格的列数
    /// @param cols 栅格列数
    void SetCols(long cols) override;
    /// @brief 获取栅格行数
    /// @return 栅格行数
    void GetRows(long *rows) override;
    /// @brief 设置栅格行数
    /// @param rows 栅格行数
    void SetRows(long rows) override;
    /// @brief 获取输出栅格的分辨率
    /// @return 输出栅格的分辨率
    void GetOutputCellSize(double_t *cellSize) override;
    /// @brief 设置输出栅格的分辨率
    /// @param outPutCellSize 输出栅格分辨率
    void SetOutputCellSize(double_t outPutCellSize) override;
  };

  // IGDALRasterReaderByPixel继承了IGDALRasterProperties和IGDALRasterReaderByPixel
  class CGDALRasterReaderByPixel : public IGDALRasterProperties, public IGDALRasterReaderByPixel {
  protected:
    long rows;
    long cols;
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
    void GetRows(long *pVal) override;
    void GetCols(long *pVal) override;
    void GetExtent(OGREnvelope **pVal) override;
    void GetCellSize(double_t *pVal) override;
    void GetBandCount(long *pVal) override;
    void GetNoData(double_t *pVal) override;
    void SetNoData(double_t pVal) override;
    void GetMinMax(bool bApproxOK, double_t *min, double_t *max, bool *pVal) override;
    void GetSpatialReference(OGRSpatialReference **pVal) override;
    void GetDataType(RasterDataType *pVal) override;
    void GetPixelValue(long col, long row, float_t *data) override;
    void ComputeRasterMinMax(bool bApproxOK, double_t *min, double_t *max, bool *pVal) override;
    void GetStatistics(bool bApproxOK, double_t *min, double_t *max, double_t *mean,
                       double_t *stddev, bool *pVal) override;
    void ComputeStatistics(bool bApproxOK, double_t *min, double_t *max, double_t *mean,
                           double_t *stddev, bool *pVal) override;
    void GetMetaData(std::vector<std::string> **pVal) override;
    void GetColorTable(std::vector<GDALColorTable> **pVal) override;
  };

  class CGDALRasterReaderByFileArray : public IGDALRasterProperties,
                                       public IGDALRasterReaderByFileArray {
  public:
    ~CGDALRasterReaderByFileArray() = default;

  protected:
    long rows;
    long cols;
    double XMin, YMin, XMax, YMax;
    double_t CellSize;
    GDALDataset *poDataset;
    GDALDataset *poSubDataset;
    GDALRasterBand *poBand;
    PBand_T pBand;
    std::vector<std::string> metadatas;
    std::string lpszPathName;
    IFileFloatArray *pData;
    long BufferSize;
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

    IFileFloatArray *GetDataBlock(long x1, long y1, long x2, long y2, long buffx, long buffy);
    bool GetDataBlock(OGRPoint LeftTop, float CellSize, int Width, int Height, float NoData);
    bool ReadDataBlock(long x1, long y1, long x2, long y2, long buffx, long buffy,
                       IFileFloatArray *pData);
    bool GetInterpolatedDataBlock(OGRPoint LeftTop, float CellSize, int Width, int Height,
                                  float NoData);
    bool GetDataBlockWithProj(OGRPoint LeftTop, float CellSize, int Width, int Height,
                              float NoData);
    bool GetInterpolatedDataBlockWithProj(OGRPoint LeftTop, float CellSize, int Width, int Height,
                                          float NoData);

  public:
    void OpenRaster(std::string lpszPathName, bool *pVal) override;
    void SetRasterBand(PBand_T band, bool *pVal) override;
    void GetPathName(std::string *pVal) override;
    void GetCurrentBand(PBand_T *pVal) override;
    void GetRows(long *pVal) override;
    void GetCols(long *pVal) override;
    void GetExtent(OGREnvelope **pVal) override;
    void GetCellSize(double_t *pVal) override;
    void GetBandCount(long *pVal) override;
    void GetNoData(double_t *pVal) override;
    void SetNoData(double_t pVal) override;
    void GetMinMax(bool bApproxOK, double_t *min, double_t *max, bool *pVal) override;
    void GetStatistics(bool bApproxOK, double_t *min, double_t *max, double_t *mean,
                       double_t *stddev, bool *pVal) override;
    void ComputeRasterMinMax(bool bApproxOK, double_t *min, double_t *max, bool *pVal) override;
    void ComputeStatistics(bool bApproxOK, double_t *min, double_t *max, double_t *mean,
                           double_t *stddev, bool *pVal) override;
    void GetSpatialReference(OGRSpatialReference **pVal) override;
    void GetMetaData(std::vector<std::string> **pVal) override;
    void GetDataType(RasterDataType *pVal) override;
    void GetColorTable(std::vector<GDALColorTable> **pVal) override;

    void GetBlockData(long x1, long y1, long x2, long y2, long buffx, long buffy,
                      IFileFloatArray **pVal) override;
    void GetInterpolatedBlockData(long x1, long y1, long x2, long y2, long BuffX, long BuffY,
                                  IFileFloatArray **pVal) override;
    void SetRefTargetSpatialReference(OGRSpatialReference *newVal) override;
    void GetBlockDataByCoord(OGRPoint *LeftTop, double_t CellSize, long Width, long Height,
                             float NoData, IFileFloatArray **pVal) override;
    void GetInterpolatedBlockDataByCoord(OGRPoint *LeftTop, double_t CellSize, long Width,
                                         long Height, float_t NoData,
                                         IFileFloatArray **pVal) override;
  };

}  // namespace spread

#endif  // !SPREADALG_INCLUDE_SPREAD_