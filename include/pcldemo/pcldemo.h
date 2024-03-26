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
    /**
     * @brief Creates a new PCLDemo
     * @param name the name to greet
     */
    PCLDemo(std::string name);

    /**
     * @brief indicate the version of program
     * @return a string containing the version info
     */
    std::string program_version() const;
  };

}  // namespace pcldemo
