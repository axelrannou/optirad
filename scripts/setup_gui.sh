#!/bin/bash

echo "=== OptiRad GUI Setup ==="
echo ""

# Detect OS
if [ -f /etc/os-release ]; then
    . /etc/os-release
    OS=$ID
else
    echo "Cannot detect OS"
    exit 1
fi

# Install OpenGL and GLFW
echo "Installing OpenGL, GLFW, and GLEW..."
if [ "$OS" = "ubuntu" ] || [ "$OS" = "debian" ]; then
    sudo apt-get update
    sudo apt-get install -y libgl1-mesa-dev libglu1-mesa-dev libglfw3-dev libglew-dev libgtest-dev libtbb-dev 
elif [ "$OS" = "fedora" ] || [ "$OS" = "rhel" ] || [ "$OS" = "centos" ]; then
    sudo dnf install -y mesa-libGL-devel mesa-libGLU-devel glfw-devel glew-devel gtest-devel tbb-devel
else
    echo "Unsupported OS: $OS"
    echo "Please manually install: OpenGL development libraries, GLFW3, and GLEW"
    exit 1
fi

echo ""
echo "Ensuring Dear ImGui docking branch..."

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMGUI_DIR="${ROOT_DIR}/external/imgui"
IMGUI_REPO="https://github.com/ocornut/imgui.git"

if [ -d "${IMGUI_DIR}/.git" ]; then
    echo "ImGui repo found at ${IMGUI_DIR}"
    git -C "${IMGUI_DIR}" fetch origin --tags
    git -C "${IMGUI_DIR}" checkout docking
    git -C "${IMGUI_DIR}" pull --ff-only origin docking
else
    echo "Cloning ImGui docking branch into ${IMGUI_DIR}"
    rm -rf "${IMGUI_DIR}"
    git clone --branch docking --single-branch "${IMGUI_REPO}" "${IMGUI_DIR}"
fi

echo "ImGui is set to branch: $(git -C "${IMGUI_DIR}" rev-parse --abbrev-ref HEAD)"