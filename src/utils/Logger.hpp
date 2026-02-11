#pragma once

#include <string>
#include <iostream>
#include <mutex>

namespace optirad {

class Logger {
private:
    static std::mutex s_mutex;
    
public:
    static void init() { /* placeholder */ }
    
    static void info(const std::string& msg) {
        std::lock_guard<std::mutex> lock(s_mutex);
        std::cout << "[INFO] " << msg << "\n";
    }
    
    static void warn(const std::string& msg) {
        std::lock_guard<std::mutex> lock(s_mutex);
        std::cout << "[WARN] " << msg << "\n";
    }
    
    static void error(const std::string& msg) {
        std::lock_guard<std::mutex> lock(s_mutex);
        std::cerr << "[ERROR] " << msg << "\n";
    }
    
    static void debug(const std::string& msg) {
        #ifdef DEBUG
        std::lock_guard<std::mutex> lock(s_mutex);
        std::cout << "[DEBUG] " << msg << "\n";
        #endif
    }
};

} // namespace optirad
