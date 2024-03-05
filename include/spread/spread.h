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

  /**
   * @brief 场强分析算法的基类
   */
  class CSpreadAnalyse {
  public:
    CSpreadAnalyse();
    virtual ~CSpreadAnalyse();

  public:
    /**
     * @brief 初始化分析环境
     */
    bool InitEnvironment();

  public:
    /**
     * @brief 分析环境参数结构体，封装了对环境参数的get与set
     */
    struct AnalyseEnvironmentST* pEnvi;

    /**
     * @brief 高程栅格路径
     */
    std::string elevationPath;

    // 封装过的根据条件获取栅格数据的对象
    // IGDALRasterReaderByFileArray*pElevs;

    /**
     * @brief 错误信息
     */
    std::string errorInfo;

    // 栅格数据对象
    // IFileFloatArray*pElevData;

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
     * @attention 默认值-9999.99，代表栅格中没有空值
     */
    double_t noData = -9999.99;

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
  protected:
    // IStations*pStations;
    float_t hm;
    std::string outputPath;
    // IFileFloatArray*pData;
    float_t offsetDB;
    bool needComputeAll;
    std::vector<double_t> crsv;
    Station stationInfo;
    double SubExtent[4];

  public:
    CFieldStrengthAnalyse();
    virtual ~CFieldStrengthAnalyse();

  protected:
    /* 子类必须实现的方法 */
    /// @brief 发射站覆盖分析
    /// @param stationInfo 台站信息
    /// @param rsv 动态数组，用于保存每个子类保存的数据，每次使用前需要留意数组中已经存在的数据
    /// @param point 待分析的点坐标
    /// @return 返回点对点之间场强分析的结果
    virtual float GetRadiuValue(Station& stationInfo, std::vector<double>& rsv, OGRPoint& point)
        = 0;

    /// @brief
    /// @param stationInfo 台站信息
    /// @param rsv
    /// @param point 待分析的点坐标
    /// @return
    virtual float GetRadiuValueRev(Station& stationInfo, std::vector<double>& rsv, OGRPoint& point)
        = 0;

    /// @brief
    /// @param stationInfo
    /// @param rsv 是一个数组，该数组里面装了台站的初始的发射功率以及其他相关信息
    virtual void PrepareReservedValues(Station& stationInfo, std::vector<double>& rsv) = 0;

    /// @brief 分析所需要的其他栅格数据
    /// @return 返回准备的状态
    virtual bool PrepareOtherData() = 0;

    /// @brief 获取模型的名称
    /// @return 返回模型的名称
    virtual std::string GetModelName() = 0;
  };

  /**
   * @brief
   * 组合分析类，是需要多种数据源以及绕射模型计算的类的父类，子类中明确了特定的绕射模型，该类中根据子类传入的rsv数组中的绕射
   * 模型调用对应的实现计算处绕射衰减
   * @attention 该类的子类模型均有绕射衰减的计算，但是绕射衰减的计算是可选的，默认值为0
   */
  class CCombineAnalyse : public CFieldStrengthAnalyse {
  protected:
    float GetRadiuValue(Station& para, std::vector<double>& rsv, double X, double Y, float Z);
    float GetRadiuValueRev(Station& para, std::vector<double>& rsv, double X, double Y, float Z);

    /// @brief
    /// 子类分析算法需要的其他数据，具体数据由该子类管理（如算法需要的特定图层，点等对象），但是给父类返回准备的结果
    /// @return 返回准备的状态
    bool PrepareOtherData();

  protected:
    // std::vector<IGDALRasterReaderByFileArray*>otherReaders;
    // std::vector<IFileFloatArray*>otherDatas;
    // std::vector<IDBLossElement*,IDBLossElement*>els;
    std::vector<std::string> otherNames;
    std::vector<std::string> otherPaths;

  protected:
    float_t GetRadiuValue(Station& stationInfo, std::vector<double_t>& rsv, OGRPoint* point);
  };

  /// @brief 自由空间传播模型，暂只实现必要方法
  class CFreeSpaceAnalyse : public CCombineAnalyse {
  protected:
    float_t GetRadiuValue(Station& stationInfo, std::vector<double_t>& rsv, OGRPoint* point);
  };

  /**
   * @brief 分析环境的参数结构体，专门管理分析环境参数
   */
  struct AnalyseEnvironmentST {
  private:
    OGRPoint* leftTop;
    double_t cellSize;
    long cols;
    long rows;
    double_t outPutCellSize;

  public:
    /// @brief 获取栅格左上角坐标
    /// @return OGRPoint*
    virtual OGRPoint* GetLeftTop() = 0;
    /// @brief 设置栅格左上角坐标
    /// @param point OGR坐标对象
    virtual void SetLetTop(OGRPoint point) = 0;

    /// @brief 获取栅格分辨率
    /// @return 栅格分辨率
    virtual double_t GetCellSize() = 0;
    /// @brief 设置栅格分辨率
    /// @param cellSize 栅格分辨率
    virtual void SetCellSize(double_t cellSize) = 0;

    /// @brief 获取栅格的列数
    /// @return 栅格列总数
    virtual long GetCols() = 0;
    /// @brief 设置栅格的列数
    /// @param cols 栅格列数
    virtual void SetCols(long cols) = 0;

    /// @brief 获取栅格行数
    /// @return 栅格行数
    virtual long GetRows() = 0;
    /// @brief 设置栅格行数
    /// @param rows 栅格行数
    virtual void SetRows(long rows) = 0;

    /// @brief 获取输出栅格的分辨率
    /// @return 输出栅格的分辨率
    virtual double_t GetOutputCellSize() = 0;
    /// @brief 设置输出栅格的分辨率
    /// @param outPutCellSize 输出栅格分辨率
    virtual void SetOutputCellSize(double_t outPutCellSize) = 0;
  };

}  // namespace spread

#endif  // !SPREADALG_INCLUDE_SPREAD_