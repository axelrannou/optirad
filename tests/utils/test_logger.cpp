#include <gtest/gtest.h>
#include "utils/Logger.hpp"
#include <sstream>

namespace optirad::tests {

class LoggerTest : public ::testing::Test {
protected:
    // Test logger basic functionality
};

TEST(LoggerTest, LogInfoDoesNotThrow) {
    EXPECT_NO_THROW({
        Logger::info("Test info message");
    });
}

TEST(LoggerTest, LogWarningDoesNotThrow) {
    EXPECT_NO_THROW({
        Logger::warn("Test warning message");
    });
}

TEST(LoggerTest, LogErrorDoesNotThrow) {
    EXPECT_NO_THROW({
        Logger::error("Test error message");
    });
}

TEST(LoggerTest, LogDebugDoesNotThrow) {
    EXPECT_NO_THROW({
        Logger::debug("Test debug message");
    });
}

TEST(LoggerTest, ConcurrentLogging) {
    // Test that logging is thread-safe
    #pragma omp parallel for num_threads(4)
    for (int i = 0; i < 100; ++i) {
        Logger::info("Thread-safe test message " + std::to_string(i));
    }
    SUCCEED();
}

TEST(LoggerTest, EmptyMessageHandled) {
    EXPECT_NO_THROW({
        Logger::info("");
        Logger::warn("");
        Logger::error("");
    });
}

TEST(LoggerTest, LongMessageHandled) {
    std::string longMsg(1000, 'x');
    EXPECT_NO_THROW({
        Logger::info(longMsg);
    });
}

} // namespace optirad::tests

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
