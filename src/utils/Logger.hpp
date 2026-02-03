#pragma once

#include <string>
#include <iostream>

namespace optirad {

class Logger {
public:
    static void init() { /* placeholder */ }
    
    static void info(const std::string& msg) {
        std::cout << "[INFO] " << msg << "\n";
    }
    
    static void warn(const std::string& msg) {
        std::cout << "[WARN] " << msg << "\n";
    }
    
    static void error(const std::string& msg) {
        std::cerr << "[ERROR] " << msg << "\n";
    }
    
    static void debug(const std::string& msg) {
        #ifdef DEBUG
        std::cout << "[DEBUG] " << msg << "\n";
        #endif
    }
};

} // namespace optirad
