#pragma once

#include <memory>
#include <string>
#include <array>
#include <vector>
#include <map>

namespace optirad {

class StructureSet;
class Structure;
struct Contour;

class RTStructParser {
public:
    RTStructParser() = default;
    ~RTStructParser() = default;
    
    // Parse RT-STRUCT file and return StructureSet
    std::unique_ptr<StructureSet> parse(const std::string& filePath);
    
private:
    // Extract ROI names from StructureSetROISequence
    std::map<int, std::string> extractROINames(void* dataset);
    
    // Parse single structure from ROIContourSequence item
    std::unique_ptr<Structure> parseStructure(void* roiContourItem, const std::map<int, std::string>& roiNames, size_t index);
    
    // Extract color from ROI Display Color tag
    bool extractColor(void* roiContourItem, std::array<uint8_t, 3>& color);
    
    // Determine structure type from name
    std::string determineType(const std::string& name);
    
    // Parse contours from ContourSequence
    std::vector<Contour> parseContours(void* contourSequence);
    
    // Parse single contour from contour item
    bool parseContour(void* contourItem, Contour& contour);
    
    // Default color palette
    static const std::vector<std::array<uint8_t, 3>> s_defaultColors;
};

}
