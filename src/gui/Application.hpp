#pragma once

#include "Window.hpp"
#include "Renderer.hpp"
#include "panels/IPanel.hpp"
#include <vector>
#include <memory>

namespace optirad {

class Application {
public:
    Application();
    ~Application();

    bool init();
    void run();
    void shutdown();

private:
    void update();
    void render();

    std::unique_ptr<Window> m_window;
    std::unique_ptr<Renderer> m_renderer;
    std::vector<std::unique_ptr<IPanel>> m_panels;
    bool m_running = false;
};

} // namespace optirad
