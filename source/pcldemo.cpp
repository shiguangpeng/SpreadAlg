#include <fmt/format.h>
#include <pcldemo/pcldemo.h>

using namespace pcldemo;

PCLDemo::PCLDemo(std::string _name) : name(std::move(_name)) {}

std::string PCLDemo::program_version() const { return fmt::format("{}", name); }
