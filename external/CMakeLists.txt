# External dependencies
include(FetchContent)

# nlohmann/json for JSON parsing (machine files, config, etc.)
FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
)
set(JSON_BuildTests OFF CACHE INTERNAL "")
FetchContent_MakeAvailable(nlohmann_json)

# Example with FetchContent:
# FetchContent_Declare(
#     eigen
#     GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
#     GIT_TAG 3.4.0
# )
# FetchContent_MakeAvailable(eigen)
#
# FetchContent_Declare(
#     spdlog
#     GIT_REPOSITORY https://github.com/gabime/spdlog.git
#     GIT_TAG v1.12.0
# )
# FetchContent_MakeAvailable(spdlog)
