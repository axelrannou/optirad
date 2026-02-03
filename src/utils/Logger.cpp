#include "Logger.hpp"
#include <iostream>

namespace optirad {

void Logger::init() {
    // TODO: Initialize spdlog
}

void Logger::log(Level level, const std::string& message) {
    std::cout << message << std::endl;
}

void Logger::debug(const std::string& message) { log(Level::Debug, message); }
void Logger::info(const std::string& message) { log(Level::Info, message); }
void Logger::warning(const std::string& message) { log(Level::Warning, message); }
void Logger::error(const std::string& message) { log(Level::Error, message); }

} // namespace optirad
