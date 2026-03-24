#pragma once

#include <string>

namespace optirad {

class IView {
public:
    virtual ~IView() = default;

    virtual std::string getName() const = 0;
    virtual void update() = 0;
    virtual void render() = 0;
    virtual void resize(int width, int height) = 0;

    void setVisible(bool visible) { m_visible = visible; }
    bool isVisible() const { return m_visible; }

protected:
    bool m_visible = true;
};

} // namespace optirad
