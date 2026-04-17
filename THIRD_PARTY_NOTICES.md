# Third-Party Notices

This project includes or is derived from third-party software components. Their respective licenses and attribution requirements remain applicable.

## matRad

Portions of OptiRad are derived from or closely based on matRad algorithms, data structures, and conventions.

- Project: matRad
- Upstream: https://github.com/e0404/matRad
- Copyright: the matRad development team within the Radiotherapy Optimization Group, German Cancer Research Center (DKFZ)
- License: BSD 3-Clause License

The matRad project states that its files are subject to the license terms in the upstream LICENSE file:
https://github.com/e0404/matRad/blob/master/LICENSE.md

OptiRad includes documentation of the specific derived and inspired components in docs/acknowledgments.md.

## Dear ImGui

OptiRad vendors Dear ImGui in external/imgui/.

- Project: Dear ImGui
- Upstream: https://github.com/ocornut/imgui
- Copyright: Omar Cornut and contributors
- License: MIT License

The full license text for the vendored copy is available in external/imgui/LICENSE.txt.

## nlohmann/json

OptiRad uses nlohmann/json via CMake FetchContent.

- Project: JSON for Modern C++
- Upstream: https://github.com/nlohmann/json
- Copyright: Niels Lohmann and contributors
- License: MIT License

## GoogleTest

OptiRad uses GoogleTest for unit tests when tests are enabled.

- Project: GoogleTest
- Upstream: https://github.com/google/googletest
- Copyright: Google and contributors
- License: BSD 3-Clause License

## DCMTK

OptiRad optionally uses DCMTK for DICOM import and export.

- Project: DCMTK
- Upstream: https://dicom.offis.de/dcmtk.php.en
- Copyright: OFFIS and contributors
- License: BSD-style license

## oneTBB

OptiRad optionally uses oneTBB for parallel RT-STRUCT parsing.

- Project: oneTBB
- Upstream: https://github.com/oneapi-src/oneTBB
- Copyright: Intel and contributors
- License: Apache License 2.0

## OpenMP

OptiRad can use OpenMP for parallel computation when available through the compiler toolchain.

- Project: OpenMP
- Upstream: https://www.openmp.org/
- License: implementation-dependent runtime licensing

## GLFW

OptiRad uses GLFW for windowing and input in the GUI build.

- Project: GLFW
- Upstream: https://www.glfw.org/
- Copyright: GLFW authors and contributors
- License: zlib/libpng License

## GLEW

OptiRad uses GLEW for OpenGL extension loading in the GUI build.

- Project: GLEW
- Upstream: http://glew.sourceforge.net/
- Copyright: GLEW authors and contributors
- License: BSD 3-Clause License / MIT-style dual licensing

## GLM

OptiRad uses GLM for mathematics in the GUI build.

- Project: OpenGL Mathematics (GLM)
- Upstream: https://github.com/g-truc/glm
- Copyright: G-Truc Creation and contributors
- License: Happy Bunny License / MIT-like permissive license

## Notes

- This file provides high-level attribution and license identification only.
- Users distributing binaries or source should also preserve the original license files and notices shipped by upstream dependencies when required.
- OptiRad itself is licensed under GPL-3.0. See the root LICENSE file.
