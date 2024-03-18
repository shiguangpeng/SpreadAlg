/**
 * all copyright reserved shigp.
 */
#ifndef INCLUDE_SPATDATA_DATATYPE_H_
#define INCLUDE_SPATDATA_DATATYPE_H_

#include <string>

using std::string;

namespace datatype {

/// @brief 描述站的结构体
typedef struct StationST {
  /// @brief 台站id
  int64_t id;
  /// @brief 东西方向坐标
  double x;
  /// @brief 南北方向坐标
  double y;
  /// @brief 经度
  double lontitude;
  /// @brief 纬度
  double latitude;
  /// @brief 台站的灵敏度
  double freqThreshold;
  /// @brief 站的高度
  double stationHeight;
  /// @brief 频率
  double frequency;
  /// @brief 功率
  double power;
  /// @brief 站所处高程
  /// @attention
  /// 我国的高程基准是青岛港验潮站的长期观测资料推算出的黄海平均海面（72.260m作为零高程面）
  double dem;
  /// @brief 站的名字
  string name;
  /// @brief 站的唯一标识
  string guid;
  /// @brief 方位角，水平方向
  double azimuth;
  /// @brief 高度角，竖直方向
  double altitude;
  /// @brief 天线直径
  double antennaDiameter;
} Station;

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
}  // namespace datatype

#endif  // INCLUDE_SPATDATA_DATATYPE_H_
