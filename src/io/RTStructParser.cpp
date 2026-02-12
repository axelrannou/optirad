#include "RTStructParser.hpp"
#include "geometry/StructureSet.hpp"
#include "geometry/Structure.hpp"
#include "utils/Logger.hpp"
#include <algorithm>

#ifdef OPTIRAD_HAS_DCMTK
#include <dcmtk/dcmdata/dctk.h>
#endif

namespace optirad {

const std::vector<std::array<uint8_t, 3>> RTStructParser::s_defaultColors = {
    {255, 0, 0},    {0, 255, 0},    {0, 0, 255},    {255, 255, 0},
    {255, 0, 255},  {0, 255, 255},  {255, 128, 0},  {128, 0, 255},
    {255, 192, 203},{0, 128, 0},
};

std::unique_ptr<StructureSet> RTStructParser::parse(const std::string& filePath) {
#ifdef OPTIRAD_HAS_DCMTK
    DcmFileFormat fileFormat;
    if (!fileFormat.loadFile(filePath.c_str()).good()) {
        Logger::error("Failed to load RT Structure Set: " + filePath);
        return nullptr;
    }
    
    DcmDataset* dataset = fileFormat.getDataset();
    
    // Extract ROI names and metadata
    auto roiNames = extractROINames(dataset);
    if (roiNames.empty()) {
        Logger::warn("No ROI names found in RT-STRUCT");
        return std::make_unique<StructureSet>();
    }
    
    // Parse ROI contours
    DcmSequenceOfItems* contourSequence = nullptr;
    if (!dataset->findAndGetSequence(DCM_ROIContourSequence, contourSequence).good()) {
        Logger::warn("No ROIContourSequence found");
        return std::make_unique<StructureSet>();
    }
    
    Logger::info("Parsing " + std::to_string(contourSequence->card()) + " structures from RT-STRUCT");
    
    auto structureSet = std::make_unique<StructureSet>();
    
    for (unsigned long i = 0; i < contourSequence->card(); ++i) {
        DcmItem* roiContourItem = contourSequence->getItem(i);
        if (!roiContourItem) continue;
        
        auto structure = parseStructure(roiContourItem, roiNames, i);
        if (structure) {
            structureSet->addStructure(std::move(structure));
        }
    }
    
    Logger::info("Loaded " + std::to_string(structureSet->getCount()) + " structures");
    
    return structureSet;
#else
    Logger::warn("DCMTK not available");
    return nullptr;
#endif
}

std::map<int, std::string> RTStructParser::extractROINames(void* datasetPtr) {
#ifdef OPTIRAD_HAS_DCMTK
    DcmDataset* dataset = static_cast<DcmDataset*>(datasetPtr);
    std::map<int, std::string> roiNames;
    
    DcmSequenceOfItems* roiSequence = nullptr;
    if (dataset->findAndGetSequence(DCM_StructureSetROISequence, roiSequence).good() && roiSequence) {
        for (unsigned long i = 0; i < roiSequence->card(); ++i) {
            DcmItem* item = roiSequence->getItem(i);
            if (!item) {
                Logger::warn("RTStructParser: null item at index " + std::to_string(i));
                continue;
            }
            
            Sint32 roiNumber = 0;
            OFString roiName;
            item->findAndGetSint32(DCM_ROINumber, roiNumber);
            item->findAndGetOFString(DCM_ROIName, roiName);
            roiNames[roiNumber] = roiName.c_str();
        }
    }
    
    return roiNames;
#else
    return {};
#endif
}

std::unique_ptr<Structure> RTStructParser::parseStructure(void* roiContourItemPtr, 
                                                           const std::map<int, std::string>& roiNames, 
                                                           size_t index) {
#ifdef OPTIRAD_HAS_DCMTK
    DcmItem* roiContourItem = static_cast<DcmItem*>(roiContourItemPtr);
    
    Sint32 refROINumber = 0;
    roiContourItem->findAndGetSint32(DCM_ReferencedROINumber, refROINumber);
    
    auto structure = std::make_unique<Structure>();
    structure->setROINumber(refROINumber);
    structure->setName(roiNames.count(refROINumber) ? roiNames.at(refROINumber) : "Unknown");
    
    // Extract color
    std::array<uint8_t, 3> color;
    if (!extractColor(roiContourItem, color)) {
        color = s_defaultColors[index % s_defaultColors.size()];
    }
    structure->setColor(color[0], color[1], color[2]);
    
    // Determine type
    structure->setType(determineType(structure->getName()));
    structure->setPriority(structure->getType() == "TARGET" ? 1 : 
                           structure->getType() == "EXTERNAL" ? 5 : 3);
    
    // Parse contours
    DcmSequenceOfItems* contourSeq = nullptr;
    if (roiContourItem->findAndGetSequence(DCM_ContourSequence, contourSeq).good()) {
        auto contours = parseContours(contourSeq);
        for (auto& contour : contours) {
            structure->addContour(contour);
        }
    }
    
    Logger::info("  - " + structure->getName() + 
                " (" + structure->getType() + 
                ", " + std::to_string(structure->getContourCount()) + " contours, " +
                "RGB: " + std::to_string(color[0]) + "/" + 
                std::to_string(color[1]) + "/" + std::to_string(color[2]) + ")");
    
    return structure;
#else
    return nullptr;
#endif
}

bool RTStructParser::extractColor(void* roiContourItemPtr, std::array<uint8_t, 3>& color) {
#ifdef OPTIRAD_HAS_DCMTK
    DcmItem* roiContourItem = static_cast<DcmItem*>(roiContourItemPtr);
    
    // Method 1: Uint16 array
    const Uint16* colorData = nullptr;
    unsigned long colorCount = 0;
    if (roiContourItem->findAndGetUint16Array(DCM_ROIDisplayColor, colorData, &colorCount).good() && colorCount >= 3) {
        color[0] = static_cast<uint8_t>(std::min(colorData[0], (Uint16)255));
        color[1] = static_cast<uint8_t>(std::min(colorData[1], (Uint16)255));
        color[2] = static_cast<uint8_t>(std::min(colorData[2], (Uint16)255));
        return true;
    }
    
    // Method 2: String parsing
    OFString colorStr0, colorStr1, colorStr2;
    if (roiContourItem->findAndGetOFString(DCM_ROIDisplayColor, colorStr0, 0).good() &&
        roiContourItem->findAndGetOFString(DCM_ROIDisplayColor, colorStr1, 1).good() &&
        roiContourItem->findAndGetOFString(DCM_ROIDisplayColor, colorStr2, 2).good()) {
        try {
            int r = std::stoi(colorStr0.c_str());
            int g = std::stoi(colorStr1.c_str());
            int b = std::stoi(colorStr2.c_str());
            color[0] = static_cast<uint8_t>(std::clamp(r, 0, 255));
            color[1] = static_cast<uint8_t>(std::clamp(g, 0, 255));
            color[2] = static_cast<uint8_t>(std::clamp(b, 0, 255));
            return true;
        } catch (...) {}
    }
    
    return false;
#else
    return false;
#endif
}

std::string RTStructParser::determineType(const std::string& name) {
    std::string upperName = name;
    std::transform(upperName.begin(), upperName.end(), upperName.begin(), ::toupper);
    
    // Target patterns
    if (upperName.find("PTV") != std::string::npos || 
        upperName.find("GTV") != std::string::npos ||
        upperName.find("CTV") != std::string::npos ||
        upperName.find("TARGET") != std::string::npos ||
        upperName.find("TUMOR") != std::string::npos) {
        return "TARGET";
    } 
    // External/body patterns
    else if (upperName.find("BODY") != std::string::npos ||
               upperName.find("EXTERNAL") != std::string::npos ||
               upperName.find("SKIN") != std::string::npos) {
        return "EXTERNAL";
    }
    // Default to OAR (Organ at Risk)
    return "OAR";
}

std::vector<Contour> RTStructParser::parseContours(void* contourSeqPtr) {
#ifdef OPTIRAD_HAS_DCMTK
    DcmSequenceOfItems* contourSeq = static_cast<DcmSequenceOfItems*>(contourSeqPtr);
    std::vector<Contour> contours;
    
    for (unsigned long c = 0; c < contourSeq->card(); ++c) {
        DcmItem* contourItem = contourSeq->getItem(c);
        if (!contourItem) continue;
        
        Contour contour;
        if (parseContour(contourItem, contour)) {
            contours.push_back(contour);
        }
    }
    
    return contours;
#else
    return {};
#endif
}

bool RTStructParser::parseContour(void* contourItemPtr, Contour& contour) {
#ifdef OPTIRAD_HAS_DCMTK
    DcmItem* contourItem = static_cast<DcmItem*>(contourItemPtr);
    
    // Method 1: Float64 array
    const Float64* contourDataF64 = nullptr;
    unsigned long dataCountF64 = 0;
    if (contourItem->findAndGetFloat64Array(DCM_ContourData, contourDataF64, &dataCountF64).good() 
        && contourDataF64 && dataCountF64 >= 3) {
        
        // Check for valid coordinate count (must be divisible by 3)
        if (dataCountF64 % 3 != 0) {
            Logger::warn("Contour data count (" + std::to_string(dataCountF64) + ") not divisible by 3");
            return false;
        }
        
        size_t numCoords = dataCountF64 / 3;
        contour.points.reserve(numCoords);
        
        for (size_t p = 0; p < numCoords; ++p) {
            contour.points.push_back({
                contourDataF64[p * 3 + 0],
                contourDataF64[p * 3 + 1],
                contourDataF64[p * 3 + 2]
            });
        }
        
        if (!contour.points.empty()) {
            contour.zPosition = contour.points[0][2];
            return true;
        }
    }
    
    // Method 2: String parsing
    DcmElement* contourDataElem = nullptr;
    if (contourItem->findAndGetElement(DCM_ContourData, contourDataElem).good() && contourDataElem) {
        unsigned long count = contourDataElem->getVM();
        if (count >= 3 && (count % 3 == 0)) {
            contour.points.reserve(count / 3);
            
            for (unsigned long p = 0; p < count; p += 3) {
                OFString xStr, yStr, zStr;
                if (contourDataElem->getOFString(xStr, p).good() &&
                    contourDataElem->getOFString(yStr, p + 1).good() &&
                    contourDataElem->getOFString(zStr, p + 2).good()) {
                    try {
                        contour.points.push_back({
                            std::stod(xStr.c_str()),
                            std::stod(yStr.c_str()),
                            std::stod(zStr.c_str())
                        });
                    } catch (...) {
                        break;
                    }
                }
            }
            
            if (!contour.points.empty()) {
                contour.zPosition = contour.points[0][2];
                return true;
            }
        }
    }
    
    return false;
#else
    return false;
#endif
}

}
