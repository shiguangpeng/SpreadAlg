/**
 * @copyright all copyright reserved
 * @author shigp
 */

#pragma once

#include <pcldemo/version.h>

#include <string>

namespace pcldemo {
/**
 * @brief A basic class for PCLDemo model
 *
 */
class PCLDemo {
    std::string name;

 public:
  explicit PCLDemo(std::string name);
};

}  // namespace pcldemo
