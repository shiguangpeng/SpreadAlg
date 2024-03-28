/**
 * @copyright all copyright reserved
 * @author shigp
 */

#include <fmt/format.h>
#include <pcldemo.h>

using pcldemo::PCLDemo;

PCLDemo::PCLDemo(std::string _name) : name(std::move(_name)) { fmt::format(_name + "PCLDemo"); }
