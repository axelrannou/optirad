#include "Timer.hpp"

namespace optirad {

void Timer::start() {
    m_start = std::chrono::high_resolution_clock::now();
}

void Timer::stop() {
    m_end = std::chrono::high_resolution_clock::now();
}

double Timer::elapsedMs() const {
    return std::chrono::duration<double, std::milli>(m_end - m_start).count();
}

double Timer::elapsedSeconds() const {
    return std::chrono::duration<double>(m_end - m_start).count();
}

} // namespace optirad
