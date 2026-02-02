#pragma once

namespace optirad {

class Renderer {
public:
    Renderer();
    ~Renderer();

    bool init();
    void shutdown();

    void beginFrame();
    void endFrame();

private:
    bool m_initialized = false;
};

} // namespace optirad
