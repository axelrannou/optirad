add_library(optirad_steering INTERFACE)

target_include_directories(optirad_steering INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/steering>
)
