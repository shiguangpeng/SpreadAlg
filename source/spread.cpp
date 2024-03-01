#include <fmt/format.h>
#include <spread/spread.h>

#include <iostream>

#include "spread/spatidatamanager.h"

using namespace spread;
using namespace spatidatamanager;

#pragma region "传播分析基类CSpreadAnalyse定义"
CSpreadAnalyse::CSpreadAnalyse() {}
CSpreadAnalyse::~CSpreadAnalyse() {}

bool CSpreadAnalyse::InitEnvironment() {
  // 通过打开的栅格数据，初始化类相关成员。
  const char* path = elevationPath.c_str();
  GDALDataset* ptr = nullptr;
  int flag = RasterDataSource::OpenRaster(path, ptr);
  if (flag == -1) {
    errorInfo = "栅格数据打开失败。";
    return false;
  }
  char** metadata = ptr->GetMetadata("");
  std::cout << "返回结果：" << *(*metadata + 1) << std::endl;
  // 结束前关闭数据集对象
  GDALClose(ptr);
  errorInfo = "栅格数据打开成功。";
  return true;
}
#pragma endregion
