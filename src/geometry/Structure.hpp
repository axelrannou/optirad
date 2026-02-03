#pragma once

#include <string>
#include <vector>
#include <cstddef>

namespace optirad {

class Structure {
public:
    void setName(const std::string& name);
    void setColor(float r, float g, float b);
    void setType(const std::string& type); // PTV, OAR, etc.

    const std::string& getName() const;
    const std::string& getType() const;

    // Mask for voxels inside structure
    void setMask(const std::vector<bool>& mask);
    const std::vector<bool>& getMask() const;
    std::vector<size_t> getVoxelIndices() const;

private:
    std::string m_name;
    std::string m_type;
    float m_color[3] = {1.0f, 0.0f, 0.0f};
    std::vector<bool> m_mask;
};

} // namespace optirad
