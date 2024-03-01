#include <fmt/format.h>
#include <spread/spread.h>

#include "gdal.h"

using namespace spread;

// ********************** CSpreadAnalyse基类定义 ********************** //

bool CSpreadAnalyse::InitEnvironment(){
    return false;
}

// Spread::Spread(std::string _name) : name(std::move(_name)) {}

// std::string Spread::program_version() const { return fmt::format("{}", name); }
// const char* Spread::test_gdal() { return GDALVersionInfo(""); }



