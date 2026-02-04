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
echo "Installing OpenGL and GLFW..."
if [ "$OS" = "ubuntu" ] || [ "$OS" = "debian" ]; then
    sudo apt-get update
    sudo apt-get install -y libgl1-mesa-dev libglu1-mesa-dev libglfw3-dev
elif [ "$OS" = "fedora" ] || [ "$OS" = "rhel" ] || [ "$OS" = "centos" ]; then
    sudo dnf install -y mesa-libGL-devel mesa-libGLU-devel glfw-devel
else
    echo "Unsupported OS: $OS"
    echo "Please manually install: OpenGL development libraries and GLFW3"
    exit 1
fi

# Download Dear ImGui
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
EXTERNAL_DIR="$PROJECT_ROOT/external"

mkdir -p "$EXTERNAL_DIR"
cd "$EXTERNAL_DIR"

if [ -d "imgui" ]; then
    echo "Dear ImGui already exists, updating..."
    cd imgui
    git pull
else
    echo "Downloading Dear ImGui..."
    git clone https://github.com/ocornut/imgui.git --branch v1.90.1 --depth 1
fi

echo ""
echo "=== Setup Complete ==="
echo "Now run:"
echo "  cd $PROJECT_ROOT/build"
echo "  cmake .."
echo "  make -j$(nproc)"
