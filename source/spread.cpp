#include <fmt/format.h>
#include <spread/spread.h>

using namespace spread;

Spread::Spread(std::string _name) : name(std::move(_name)) {}

std::string Spread::program_version() const { return fmt::format("{}", name); }
