#pragma once

#include <chrono>
#include <string>

namespace optirad {

class Timer {
public:
    Timer() : m_start(std::chrono::high_resolution_clock::now()),
              m_end(std::chrono::high_resolution_clock::now()) {}
    
    void start();
    void stop();
    double elapsedMs() const;
    double elapsedSeconds() const;

private:
    std::chrono::high_resolution_clock::time_point m_start;
    std::chrono::high_resolution_clock::time_point m_end;
};

} // namespace optirad
