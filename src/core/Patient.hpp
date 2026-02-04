#pragma once

#include <string>

namespace optirad {

class Patient {
public:
    void setName(const std::string& name) { m_name = name; }
    void setID(const std::string& id) { m_id = id; }
    
    const std::string& getName() const { return m_name; }
    const std::string& getID() const { return m_id; }

private:
    std::string m_name;
    std::string m_id;
};

} // namespace optirad
