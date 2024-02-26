#pragma once

#include <string>
#include <spread/version.h>

namespace spread {
  /**
   * @brief A basic class for spread model 
   * 
   */
  class Spread {
    std::string name;

  public:
    /**
     * @brief Creates a new spread
     * @param name the name to greet
     */
    Spread(std::string name);

    /**
     * @brief indicate the version of program
     * @return a string containing the version info
     */
    std::string program_version() const;
  };

}  // namespace spread
