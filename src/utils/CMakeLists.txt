add_library(optirad_utils
    Logger.cpp
    Config.cpp
    Timer.cpp
    MathUtils.cpp
    Interpolation.cpp
)

target_include_directories(optirad_utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# target_link_libraries(optirad_utils PUBLIC spdlog::spdlog nlohmann_json::nlohmann_json)
