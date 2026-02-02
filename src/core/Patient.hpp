#pragma once

#include <string>

namespace optirad {

class Patient {
public:
    void setId(const std::string& id);
    void setName(const std::string& name);

    const std::string& getId() const;
    const std::string& getName() const;

private:
    std::string m_id;
    std::string m_name;
};

} // namespace optirad
