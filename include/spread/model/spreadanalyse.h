/**
 * @copyright ALL COPYRIGH RESERVED BY SHIGP
 * @author shigp
 */

#ifndef INCLUDE_SPREAD_MODEL_SPREADANALYSE_H_
#define INCLUDE_SPREAD_MODEL_SPREADANALYSE_H_

#include <spatdata/datatype.h>
#include <spatdata/spatdata.h>

#include <string>
#include <vector>

using datatype::Station;
using spatdata::IAnalyseEnvironment;
using spatdata::IFileFloatArray;
using spatdata::IGDALRasterReaderByFileArray;
using std::string;
using std::vector;

namespace spread {
namespace spreadanalyse {
struct AnalyseEnvironment : public IAnalyseEnvironment {
 public:
  AnalyseEnvironment();
  // ~AnalyseEnvironment();

 private:
  OGRPoint *leftTop;
  double_t cellSize;
  int64_t cols;
  int64_t rows;
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
  void GetCols(int64_t *cols) override;
  /// @brief 设置栅格的列数
  /// @param cols 栅格列数
  void SetCols(int64_t cols) override;
  /// @brief 获取栅格行数
  /// @return 栅格行数
  void GetRows(int64_t *rows) override;
  /// @brief 设置栅格行数
  /// @param rows 栅格行数
  void SetRows(int64_t rows) override;
  /// @brief 获取输出栅格的分辨率
  /// @return 输出栅格的分辨率
  void GetOutputCellSize(double_t *cellSize) override;
  /// @brief 设置输出栅格的分辨率
  /// @param outPutCellSize 输出栅格分辨率
  void SetOutputCellSize(double_t outPutCellSize) override;
};

// struct ISpreadAnalyse {
//   virtual bool InitEnvironment();
// };

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
  int64_t cols;

  /**
   * @brief 栅格行数
   */
  int64_t rows;

  /**
   * @brief 栅格大小
   */
  double cellSize;

  /**
   * @brief 最小经度位置（左上角）
   */
  double xMin;

  /**
   * @brief 最大纬度位置（左上角）
   */
  double yMax;

  /**
   * @brief 栅格无数据（或无有效值）的填充值
   */
  double noData;

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

}  // namespace spreadanalyse

}  // namespace spread
#endif  // INCLUDE_SPREAD_MODEL_SPREADANALYSE_H_
