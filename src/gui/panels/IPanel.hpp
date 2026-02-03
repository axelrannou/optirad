#pragma once

#include <string>

namespace optirad {

class IPanel {
public:
    virtual ~IPanel() = default;

    virtual std::string getName() const = 0;
    virtual void update() = 0;
    virtual void render() = 0;

    void setVisible(bool visible) { m_visible = visible; }
    bool isVisible() const { return m_visible; }

protected:
    bool m_visible = true;
};

} // namespace optirad
