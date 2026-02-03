#pragma once

#include <chrono>
#include <string>

namespace optirad {

class Timer {
public:
    void start();
    void stop();
    double elapsedMs() const;
    double elapsedSeconds() const;

private:
    std::chrono::high_resolution_clock::time_point m_start;
    std::chrono::high_resolution_clock::time_point m_end;
};

} // namespace optirad
