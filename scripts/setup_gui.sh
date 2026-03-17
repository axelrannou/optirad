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