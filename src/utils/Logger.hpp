#pragma once

#include <string>

namespace optirad {

class Logger {
public:
    enum class Level { Debug, Info, Warning, Error };

    static void init();
    static void log(Level level, const std::string& message);
    static void debug(const std::string& message);
    static void info(const std::string& message);
    static void warning(const std::string& message);
    static void error(const std::string& message);
};

} // namespace optirad
