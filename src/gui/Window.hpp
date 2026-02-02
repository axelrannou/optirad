#pragma once

#include <string>

namespace optirad {

class Window {
public:
    Window();
    ~Window();

    bool create(const std::string& title, int width, int height);
    void destroy();

    void pollEvents();
    void swapBuffers();
    bool shouldClose() const;

    int getWidth() const;
    int getHeight() const;

private:
    void* m_handle = nullptr; // GLFWwindow*
    int m_width = 0;
    int m_height = 0;
};

} // namespace optirad
