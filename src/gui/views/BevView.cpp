#include "BevView.hpp"
#include "../panels/DeliveryPanel.hpp"
#include "core/Aperture.hpp"
#include "core/FluenceMap.hpp"
#include "geometry/MathUtils.hpp"
#include <imgui.h>
#include <cmath>
#include <cstdio>
#include <algorithm>

namespace optirad {

BevView::BevView(GuiAppState& state) : m_state(state) {
    glGenTextures(1, &m_fluenceTexID);
}

BevView::~BevView() {
    if (m_fluenceTexID) glDeleteTextures(1, &m_fluenceTexID);
}

void BevView::jetColormap(float t, unsigned char& r, unsigned char& g, unsigned char& b) {
    t = std::clamp(t, 0.0f, 1.0f);
    float r4 = std::clamp(1.5f - std::abs(t - 0.75f) * 4.0f, 0.0f, 1.0f);
    float g4 = std::clamp(1.5f - std::abs(t - 0.50f) * 4.0f, 0.0f, 1.0f);
    float b4 = std::clamp(1.5f - std::abs(t - 0.25f) * 4.0f, 0.0f, 1.0f);
    r = static_cast<unsigned char>(r4 * 255);
    g = static_cast<unsigned char>(g4 * 255);
    b = static_cast<unsigned char>(b4 * 255);
}

ImVec2 BevView::bevToScreen(double bevX, double bevZ,
                             ImVec2 imgMin, ImVec2 imgMax) const {
    double rangeX = m_bevMaxX - m_bevMinX;
    double rangeZ = m_bevMaxZ - m_bevMinZ;
    if (rangeX <= 0) rangeX = 1;
    if (rangeZ <= 0) rangeZ = 1;

    // UV in [0,1]
    float u = static_cast<float>((bevX - m_bevMinX) / rangeX);
    float v = 1.0f - static_cast<float>((bevZ - m_bevMinZ) / rangeZ);

    // Apply centered zoom
    float halfU = 0.5f / m_zoom;
    float halfV = 0.5f / m_zoom;
    float uvMinX = 0.5f - halfU;
    float uvMaxX = 0.5f + halfU;
    float uvMinY = 0.5f - halfV;
    float uvMaxY = 0.5f + halfV;

    float sx = (u - uvMinX) / (uvMaxX - uvMinX);
    float sy = (v - uvMinY) / (uvMaxY - uvMinY);

    return ImVec2(imgMin.x + sx * (imgMax.x - imgMin.x),
                  imgMin.y + sy * (imgMax.y - imgMin.y));
}

void BevView::render() {
    if (!m_visible) return;

    ImGui::Begin("BEV", &m_visible, ImGuiWindowFlags_HorizontalScrollbar);
    m_state.delivery.synchronize(m_state);

    const auto* selected = m_state.doseStore.getSelected();
    const int selectedId = selected ? selected->id : -1;
    if (selectedId != m_selectedDoseId) {
        m_selectedDoseId = selectedId;
        m_segmentIndex = 0;
        m_lastDeliveryTime = -1.0;
        m_needsUpdate = true;
    }
    if (m_state.delivery.available() && m_state.delivery.isVmat()) {
        renderArcMotion();
        ImGui::End();
        return;
    }

    if (!m_state.stfGenerated()) {
        ImGui::TextDisabled("Generate STF first to view BEV.");
        ImGui::End();
        return;
    }

    renderControls();
    ImGui::Separator();
    renderBevContent();

    ImGui::End();
}

void BevView::renderArcMotion() {
    auto& delivery = m_state.delivery;
    DeliveryPanel::renderControls(delivery);
    const auto& current = delivery.frame();
    const double extentX = delivery.extentX(), extentZ = delivery.extentZ();
    renderSceneControls();
    if (m_showAnatomy) {
        std::vector<const Aperture*> apertures;
        if (current.hasAperture) apertures.push_back(&current.aperture);
        renderScene(current.beam, apertures, current.boundaries, extentX, extentZ);
        return;
    }
    const ImVec2 available = ImGui::GetContentRegionAvail();
    if (available.x <= 0.0f || available.y <= 0.0f) return;
    const ImVec2 origin = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("##arc_bev_zoom", available);
    handleZoom(m_zoom);
    const float scale = std::min(available.x / static_cast<float>(2.0 * extentX),
                                 available.y / static_cast<float>(2.0 * extentZ)) * m_zoom;
    auto screen = [&](double lateral, double longitudinal) {
        return ImVec2(origin.x + available.x * 0.5f + static_cast<float>(lateral) * scale,
                      origin.y + available.y * 0.5f - static_cast<float>(longitudinal) * scale);
    };
    auto* draw = ImGui::GetWindowDrawList();
    draw->PushClipRect(origin, ImVec2(origin.x + available.x, origin.y + available.y), true);
    draw->AddRectFilled(screen(-extentX, extentZ), screen(extentX, -extentZ), IM_COL32(24, 30, 34, 255));
    for (size_t pair = 0; pair < current.aperture.bankA.size(); ++pair) {
        const double left = current.aperture.bankA[pair];
        const double right = current.aperture.bankB[pair];
        const double bottom = current.boundaries[pair], top = current.boundaries[pair + 1];
        draw->AddRectFilled(screen(left, top), screen(right, bottom), IM_COL32(45, 178, 150, 210));
        for (const auto& edges : {std::pair<double, double>{-extentX, left}, {right, extentX}}) {
            draw->AddRectFilled(screen(edges.first, top), screen(edges.second, bottom), IM_COL32(100, 114, 126, 255));
            draw->AddRect(screen(edges.first, top), screen(edges.second, bottom), IM_COL32(35, 43, 48, 255));
        }
    }
    draw->PopClipRect();
}

void BevView::renderControls() {
    int numBeams = static_cast<int>(m_state.stf->getCount());
    if (numBeams == 0) return;

    DeliveryPanel::renderControls(m_state.delivery);
    if (m_state.delivery.available()) {
        const auto& frame = m_state.delivery.frame();
        if (m_lastDeliveryTime != m_state.delivery.time() || m_state.delivery.playing) {
            m_beamIndex = static_cast<int>(frame.beamIndex);
            m_segmentIndex = frame.segment;
            m_needsUpdate = true;
        }
        m_lastDeliveryTime = m_state.delivery.time();
    }

    // Clamp beam index
    if (m_beamIndex >= numBeams) m_beamIndex = numBeams - 1;
    if (m_beamIndex < 0) m_beamIndex = 0;

    // Beam slider
    bool changed = false;
    const auto* beam = m_state.stf->getBeam(static_cast<size_t>(m_beamIndex));
    std::string beamLabel = "Beam " + std::to_string(m_beamIndex) +
        " (GA: " + std::to_string(static_cast<int>(beam->getGantryAngle())) + ")";

    if (!m_state.delivery.available() &&
        ImGui::SliderInt("Beam", &m_beamIndex, 0, numBeams - 1, beamLabel.c_str())) {
        changed = true;
        beam = m_state.stf->getBeam(static_cast<size_t>(m_beamIndex));
    }

    // Segment slider (only if viewing a deliverable dose with leaf sequencing data)
    const bool isDeliverableDose = m_state.deliveredSequences() != nullptr;
    if (isDeliverableDose && m_state.leafSequenceDone() &&
        m_beamIndex < static_cast<int>(m_state.leafSequences.size())) {
        const auto& seq = m_state.leafSequences[static_cast<size_t>(m_beamIndex)];
        int numSegs = static_cast<int>(seq.segments.size());
        if (numSegs > 0) {
            m_segmentIndex = std::clamp(m_segmentIndex, -1, numSegs - 1);
            // -1 = "All", 0..numSegs-1 = individual segments
            std::string segLabel;
            if (m_segmentIndex < 0) {
                segLabel = "All (" + std::to_string(numSegs) + " segments)";
            } else {
                char muBuf[32];
                std::snprintf(muBuf, sizeof(muBuf), "%.1f",
                    seq.segments[static_cast<size_t>(m_segmentIndex)].weight);
                segLabel = "Segment " + std::to_string(m_segmentIndex) +
                    " (" + muBuf + " MU)";
            }
            if (ImGui::SliderInt("Segment", &m_segmentIndex, -1, numSegs - 1,
                                 segLabel.c_str())) {
                changed = true;
                m_state.delivery.playing = false;
                if (m_segmentIndex >= 0 && m_state.delivery.available() &&
                    seq.segments[static_cast<size_t>(m_segmentIndex)].weight > 0.0) {
                    m_state.delivery.selectSegment(static_cast<size_t>(m_beamIndex), m_segmentIndex);
                    m_lastDeliveryTime = m_state.delivery.time();
                }
            }
        }
    } else {
        m_segmentIndex = -1;
    }

    if (changed) m_needsUpdate = true;

    renderSceneControls();

    // Info line — show different details for optimization vs deliverable dose
    if (beam) {
        const bool viewingDeliverable = m_state.deliveredSequences() != nullptr;
        const auto& fs = beam->getFieldSize();
        if (viewingDeliverable && m_state.plan) {
            const auto& mlc = m_state.plan->getMachine().getGeometry();
            std::string leafInfo;
            if (mlc.leafWidths.size() == 2)
                leafInfo = std::to_string(static_cast<int>(mlc.leafWidths[0])) + "/"
                         + std::to_string(static_cast<int>(mlc.leafWidths[1])) + " mm";
            else if (!mlc.leafWidths.empty())
                leafInfo = std::to_string(static_cast<int>(mlc.leafWidths[0])) + " mm";
            else
                leafInfo = "?";
            ImGui::Text("Field: %.0f x %.0f mm | Leaf: %s | Res: %.1f mm | SAD: %.0f mm",
                        fs[0], fs[1], leafInfo.c_str(), mlc.leafPositionResolution, beam->getSAD());
            if (!mlc.mlcType.empty())
                ImGui::Text("MLC: %s", mlc.mlcType.c_str());
        } else {
            ImGui::Text("Field: %.0f x %.0f mm | Bixel: %.1f mm | SAD: %.0f mm",
                        fs[0], fs[1], beam->getBixelWidth(), beam->getSAD());
            if (m_state.plan) {
                const auto& mlcName = m_state.plan->getMachine().getGeometry().mlcType;
                if (!mlcName.empty())
                    ImGui::Text("MLC: %s", mlcName.c_str());
            }
        }
    }
}

void BevView::renderBevContent() {
    const auto* beam = m_state.stf->getBeam(static_cast<size_t>(m_beamIndex));
    if (!beam || beam->getNumOfRays() == 0) return;

    // Compute BEV extent from ray positions (with half-bixel margin)
    double bw = beam->getBixelWidth();
    double xMin = 1e30, xMax = -1e30, zMin = 1e30, zMax = -1e30;
    for (size_t i = 0; i < beam->getNumOfRays(); ++i) {
        const auto& pos = beam->getRay(i)->getRayPosBev();
        xMin = std::min(xMin, pos[0]);
        xMax = std::max(xMax, pos[0]);
        zMin = std::min(zMin, pos[2]);
        zMax = std::max(zMax, pos[2]);
    }
    m_bevMinX = xMin - bw;
    m_bevMaxX = xMax + bw;
    m_bevMinZ = zMin - bw;
    m_bevMaxZ = zMax + bw;

    if (m_showAnatomy) {
        std::vector<const Aperture*> apertures;
        std::vector<double> boundaries;
        if (m_state.deliveredSequences() && m_state.leafSequenceDone() &&
            m_beamIndex < static_cast<int>(m_state.leafSequences.size())) {
            const auto& sequence = m_state.leafSequences[static_cast<size_t>(m_beamIndex)];
            for (size_t segment = 0; segment < sequence.segments.size(); ++segment) {
                if (m_segmentIndex < 0 || segment == static_cast<size_t>(m_segmentIndex))
                    apertures.push_back(&sequence.segments[segment]);
            }
            boundaries = sequence.leafPairBoundariesZ;
            if (boundaries.empty() && !apertures.empty()) {
                for (size_t pair = 0; pair <= apertures.front()->bankA.size(); ++pair)
                    boundaries.push_back(zMin - bw * 0.5 + pair * bw);
            }
            if (apertures.size() > 1)
                ImGui::TextDisabled("All segments: aperture outlines. Select a segment to see leaves.");
        }
        double extentX = std::max(std::abs(m_bevMinX), std::abs(m_bevMaxX)) + 15.0;
        double extentZ = std::max(std::abs(m_bevMinZ), std::abs(m_bevMaxZ)) + 5.0;
        for (const auto* aperture : apertures) {
            for (double tip : aperture->bankA) extentX = std::max(extentX, std::abs(tip) + 15.0);
            for (double tip : aperture->bankB) extentX = std::max(extentX, std::abs(tip) + 15.0);
        }
        for (double boundary : boundaries) extentZ = std::max(extentZ, std::abs(boundary) + 5.0);
        const auto& frame = m_state.delivery.frame();
        if (m_state.delivery.available() && m_segmentIndex >= 0 && frame.hasAperture &&
            frame.segment == m_segmentIndex && frame.beamIndex == static_cast<size_t>(m_beamIndex)) {
            apertures = {&frame.aperture};
            boundaries = frame.boundaries;
        }
        renderScene(m_state.delivery.available() ? frame.beam : *beam,
                    apertures, boundaries, extentX, extentZ);
        return;
    }

    // Update fluence texture
    int ver = m_state.doseStore.version();
    if (ver != m_doseVersion) {
        m_needsUpdate = true;
        m_doseVersion = ver;
    }
    if (m_needsUpdate && m_state.optimizationDone()) {
        updateFluenceTexture();
        m_needsUpdate = false;
    }

    // Available region for the BEV image
    ImVec2 avail = ImGui::GetContentRegionAvail();
    if (avail.x <= 0 || avail.y <= 0) return;

    // Maintain aspect ratio
    double bevW = m_bevMaxX - m_bevMinX;
    double bevH = m_bevMaxZ - m_bevMinZ;
    if (bevW <= 0 || bevH <= 0) return;

    float aspect = static_cast<float>(bevW / bevH);
    float imgW = avail.x;
    float imgH = avail.y;
    if (imgW / imgH > aspect) {
        imgW = imgH * aspect;
    } else {
        imgH = imgW / aspect;
    }

    // Center the image in available region
    ImVec2 cursorPos = ImGui::GetCursorScreenPos();
    float offsetX = (avail.x - imgW) * 0.5f;
    float offsetY = (avail.y - imgH) * 0.5f;
    ImVec2 imgMin(cursorPos.x + offsetX, cursorPos.y + offsetY);
    ImVec2 imgMax(imgMin.x + imgW, imgMin.y + imgH);

    // Draw fluence texture (if available)
    ImDrawList* dl = ImGui::GetWindowDrawList();
    if (m_state.optimizationDone() && m_fluenceTexID && m_texWidth > 0) {
        // Determine Z range for texture placement
        double bw = beam->getBixelWidth();
        double rayXMin = m_bevMinX + bw;
        double rayXMax = m_bevMaxX - bw;
        double texZMin, texZMax;

        bool useLeafBounds = false;
        {
            if (m_state.deliveredSequences() && m_state.leafSequenceDone() &&
                m_beamIndex < static_cast<int>(m_state.leafSequences.size())) {
                const auto& seq = m_state.leafSequences[static_cast<size_t>(m_beamIndex)];
                if (!seq.leafPairBoundariesZ.empty()) {
                    texZMin = seq.leafPairBoundariesZ.front();
                    texZMax = seq.leafPairBoundariesZ.back();
                    useLeafBounds = true;
                }
            }
        }
        if (!useLeafBounds) {
            double rayZMin = m_bevMinZ + bw;
            double rayZMax = m_bevMaxZ - bw;
            texZMin = rayZMin - bw * 0.5;
            texZMax = rayZMax + bw * 0.5;
        }

        ImVec2 texTL = bevToScreen(rayXMin - bw * 0.5, texZMax, imgMin, imgMax);
        ImVec2 texBR = bevToScreen(rayXMax + bw * 0.5, texZMin, imgMin, imgMax);
        dl->PushClipRect(imgMin, imgMax, true);
        dl->AddImage(
            reinterpret_cast<void*>(static_cast<intptr_t>(m_fluenceTexID)),
            texTL, texBR, ImVec2(0, 1), ImVec2(1, 0));
        dl->PopClipRect();
    } else {
        // Dark background
        dl->AddRectFilled(imgMin, imgMax, IM_COL32(30, 30, 30, 255));
    }

    // Draw field boundary
    dl->AddRect(imgMin, imgMax, IM_COL32(200, 200, 200, 180), 0.0f, 0, 1.5f);

    // Draw MLC leaves (only for deliverable doses)
    const bool isDeliverableDose = m_state.deliveredSequences() != nullptr;
    if (isDeliverableDose && m_state.leafSequenceDone() &&
        m_beamIndex < static_cast<int>(m_state.leafSequences.size())) {
        renderMlcLeaves(dl, imgMin, imgMax);
    }

    ImGui::SetCursorScreenPos(imgMin);
    ImGui::InvisibleButton("##bev_interact", ImVec2(imgW, imgH));
    handleZoom(m_zoom);
}

void BevView::updateFluenceTexture() {
    const auto* beam = m_state.stf->getBeam(static_cast<size_t>(m_beamIndex));
    if (!beam) return;

    // Compute global bixel offset for this beam
    size_t offset = 0;
    for (int bi = 0; bi < m_beamIndex; ++bi) {
        offset += m_state.stf->getBeam(static_cast<size_t>(bi))->getTotalNumOfBixels();
    }

    // Determine if viewing deliverable dose (leaf-pair resolution) or optimization dose (bixel resolution)
    const bool isDeliverable = m_state.deliveredSequences() != nullptr;

    // For deliverable dose: use leaf-pair-resolution fluence from the sequencer result
    // For optimization dose: use bixel-resolution fluence from optimizer weights
    int texRows = 0, texCols = 0;
    double maxFluence = 1.0;
    std::vector<unsigned char> pixels;

    if (isDeliverable && m_state.leafSequenceDone() &&
        m_beamIndex < static_cast<int>(m_state.leafSequences.size())) {
        const auto& seq = m_state.leafSequences[static_cast<size_t>(m_beamIndex)];
        if (!seq.leafPairFluence.empty() && seq.leafPairFluenceCols > 0) {
            texCols = seq.leafPairFluenceCols;
            int numLP = static_cast<int>(seq.leafPairFluence.size()) / texCols;

            maxFluence = *std::max_element(seq.leafPairFluence.begin(), seq.leafPairFluence.end());
            if (maxFluence <= 0) maxFluence = 1.0;

            // Build a physically-correct texture: each pixel row spans a uniform
            // physical Z step equal to the minimum leaf width.  This ensures that
            // wider leaf pairs (e.g. 5mm HD120 outer) occupy proportionally more
            // pixel rows than narrow ones (e.g. 2.5mm inner), so the texture
            // aligns correctly with the MLC leaf boundaries drawn by renderMlcLeaves().
            const bool hasBounds = (!seq.leafPairBoundariesZ.empty() &&
                static_cast<int>(seq.leafPairBoundariesZ.size()) >= numLP + 1);

            if (hasBounds) {
                // Find minimum leaf-pair width to use as the pixel step
                double minLeafW = 1e30;
                for (int li = 0; li < numLP; ++li)
                    minLeafW = std::min(minLeafW,
                        seq.leafPairBoundariesZ[li + 1] - seq.leafPairBoundariesZ[li]);
                if (minLeafW <= 0.0 || minLeafW > 1e20) minLeafW = 2.5;

                double zBot = seq.leafPairBoundariesZ.front();
                double zTop = seq.leafPairBoundariesZ.back();
                texRows = std::max(1, static_cast<int>(std::round((zTop - zBot) / minLeafW)));
                double pixH = (zTop - zBot) / texRows;

                pixels.assign(static_cast<size_t>(texRows) * texCols * 4, 0);

                for (int pr = 0; pr < texRows; ++pr) {
                    double zCenter = zBot + (pr + 0.5) * pixH;
                    // Identify which leaf pair owns this Z position
                    int li = -1;
                    for (int lp = 0; lp < numLP; ++lp) {
                        if (zCenter >= seq.leafPairBoundariesZ[lp] &&
                            zCenter <  seq.leafPairBoundariesZ[lp + 1]) {
                            li = lp;
                            break;
                        }
                    }
                    if (li < 0) continue;

                    for (int c = 0; c < texCols; ++c) {
                        double val = seq.leafPairFluence[static_cast<size_t>(li) * texCols + c];
                        float t = static_cast<float>(val / maxFluence);
                        size_t idx = (static_cast<size_t>(pr) * texCols + c) * 4;
                        if (val > 0.001 * maxFluence) {
                            jetColormap(t, pixels[idx], pixels[idx + 1], pixels[idx + 2]);
                            pixels[idx + 3] = 200;
                        } else {
                            pixels[idx] = pixels[idx + 1] = pixels[idx + 2] = 0;
                            pixels[idx + 3] = 0;
                        }
                    }
                }
            } else {
                // Fallback (no boundary data): one row per leaf pair
                texRows = numLP;
                pixels.resize(static_cast<size_t>(texRows) * texCols * 4);
                for (int r = 0; r < texRows; ++r) {
                    for (int c = 0; c < texCols; ++c) {
                        double val = seq.leafPairFluence[static_cast<size_t>(r) * texCols + c];
                        float t = static_cast<float>(val / maxFluence);
                        size_t idx = (static_cast<size_t>(r) * texCols + c) * 4;
                        if (val > 0.001 * maxFluence) {
                            jetColormap(t, pixels[idx], pixels[idx + 1], pixels[idx + 2]);
                            pixels[idx + 3] = 200;
                        } else {
                            pixels[idx] = pixels[idx + 1] = pixels[idx + 2] = 0;
                            pixels[idx + 3] = 0;
                        }
                    }
                }
            }
        }
    }

    // Fallback: bixel-resolution texture from optimizer weights
    if (pixels.empty()) {
        const auto& weights = (isDeliverable && !m_state.deliverableWeights.empty())
            ? m_state.deliverableWeights : m_state.optimizedWeights;

        auto fm = FluenceMap::fromBeamWeights(*beam, weights, offset);
        texRows = fm.getNumRows();
        texCols = fm.getNumCols();
        if (texRows == 0 || texCols == 0) return;

        maxFluence = fm.getMaxFluence();
        if (maxFluence <= 0) maxFluence = 1.0;

        pixels.resize(static_cast<size_t>(texRows) * texCols * 4);
        for (int r = 0; r < texRows; ++r) {
            for (int c = 0; c < texCols; ++c) {
                double val = fm.getValue(r, c);
                float t = static_cast<float>(val / maxFluence);
                size_t idx = (static_cast<size_t>(r) * texCols + c) * 4;
                if (val > 0.001 * maxFluence) {
                    jetColormap(t, pixels[idx], pixels[idx + 1], pixels[idx + 2]);
                    pixels[idx + 3] = 200;
                } else {
                    pixels[idx] = pixels[idx + 1] = pixels[idx + 2] = 0;
                    pixels[idx + 3] = 0;
                }
            }
        }
    }

    if (texRows == 0 || texCols == 0) return;

    m_texWidth = texCols;
    m_texHeight = texRows;

    // Upload to OpenGL
    glBindTexture(GL_TEXTURE_2D, m_fluenceTexID);
    if (m_texWidth != m_lastTexWidth || m_texHeight != m_lastTexHeight) {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, m_texWidth, m_texHeight,
                     0, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        m_lastTexWidth = m_texWidth;
        m_lastTexHeight = m_texHeight;
    } else {
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, m_texWidth, m_texHeight,
                        GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
    }
    glBindTexture(GL_TEXTURE_2D, 0);
}

void BevView::renderMlcLeaves(ImDrawList* dl, ImVec2 imgMin, ImVec2 imgMax) {
    const auto& seq = m_state.leafSequences[static_cast<size_t>(m_beamIndex)];
    if (seq.segments.empty()) return;

    const auto* beam = m_state.stf->getBeam(static_cast<size_t>(m_beamIndex));
    if (!beam || beam->getNumOfRays() == 0) return;

    // Determine which segments to draw
    std::vector<size_t> segIndicesToDraw;
    if (m_segmentIndex < 0) {
        for (size_t i = 0; i < seq.segments.size(); ++i)
            segIndicesToDraw.push_back(i);
    } else if (m_segmentIndex < static_cast<int>(seq.segments.size())) {
        segIndicesToDraw.push_back(static_cast<size_t>(m_segmentIndex));
    }

    ImU32 leafClosedColor = IM_COL32(80, 80, 100, 160);
    ImU32 leafBorderColor = IM_COL32(140, 140, 160, 200);

    // Use physical leaf pair boundaries if available, otherwise fall back to bixel width
    const bool hasLeafBounds = !seq.leafPairBoundariesZ.empty();

    // Fallback: compute bixel-based Z origin
    double bw = beam->getBixelWidth();
    double originZ = 0.0;
    if (!hasLeafBounds) {
        originZ = 1e30;
        for (size_t i = 0; i < beam->getNumOfRays(); ++i)
            originZ = std::min(originZ, beam->getRay(i)->getRayPosBev()[2]);
    }

    for (size_t si : segIndicesToDraw) {
        const auto& frame = m_state.delivery.frame();
        const auto& seg = m_state.delivery.available() && m_segmentIndex >= 0 && frame.hasAperture &&
                          frame.segment == m_segmentIndex && frame.beamIndex == static_cast<size_t>(m_beamIndex) ?
                          frame.aperture : seq.segments[si];
        int nLeafPairs = static_cast<int>(seg.bankA.size());

        for (int r = 0; r < nLeafPairs; ++r) {
            double zBot, zTop;
            if (hasLeafBounds && r + 1 < static_cast<int>(seq.leafPairBoundariesZ.size())) {
                zBot = seq.leafPairBoundariesZ[r];
                zTop = seq.leafPairBoundariesZ[r + 1];
            } else {
                double zCenter = originZ + r * bw;
                zBot = zCenter - bw * 0.5;
                zTop = zCenter + bw * 0.5;
            }

            double leafA = seg.bankA[r];
            double leafB = seg.bankB[r];

            // Left leaf: from field left edge to bankA position
            ImVec2 p0 = bevToScreen(m_bevMinX, zBot, imgMin, imgMax);
            ImVec2 p1 = bevToScreen(leafA, zTop, imgMin, imgMax);
            ImVec2 leftMin(std::min(p0.x, p1.x), std::min(p0.y, p1.y));
            ImVec2 leftMax(std::max(p0.x, p1.x), std::max(p0.y, p1.y));

            leftMin.x = std::max(leftMin.x, imgMin.x);
            leftMin.y = std::max(leftMin.y, imgMin.y);
            leftMax.x = std::min(leftMax.x, imgMax.x);
            leftMax.y = std::min(leftMax.y, imgMax.y);

            if (leftMax.x > leftMin.x && leftMax.y > leftMin.y) {
                dl->AddRectFilled(leftMin, leftMax, leafClosedColor);
                dl->AddRect(leftMin, leftMax, leafBorderColor, 0.0f, 0, 0.5f);
            }

            // Right leaf: from bankB position to field right edge
            p0 = bevToScreen(leafB, zBot, imgMin, imgMax);
            p1 = bevToScreen(m_bevMaxX, zTop, imgMin, imgMax);
            ImVec2 rightMin(std::min(p0.x, p1.x), std::min(p0.y, p1.y));
            ImVec2 rightMax(std::max(p0.x, p1.x), std::max(p0.y, p1.y));

            rightMin.x = std::max(rightMin.x, imgMin.x);
            rightMin.y = std::max(rightMin.y, imgMin.y);
            rightMax.x = std::min(rightMax.x, imgMax.x);
            rightMax.y = std::min(rightMax.y, imgMax.y);

            if (rightMax.x > rightMin.x && rightMax.y > rightMin.y) {
                dl->AddRectFilled(rightMin, rightMax, leafClosedColor);
                dl->AddRect(rightMin, rightMax, leafBorderColor, 0.0f, 0, 0.5f);
            }
        }
    }
}

void BevView::renderSceneControls() {
    ImGui::Checkbox("3D anatomy", &m_showAnatomy);
    ImGui::SameLine();
    if (ImGui::Button("Reset zoom")) {
        m_sceneZoom = 1.0f;
        m_zoom = 1.0f;
    }
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Reset magnification. Scroll over the view to zoom; the view stays source-aligned.");
}

void BevView::handleZoom(float& zoom) {
    if (ImGui::IsItemHovered()) {
        zoom *= std::exp(ImGui::GetIO().MouseWheel * 0.1f);
        zoom = std::clamp(zoom, 0.25f, 10.0f);
    }
}

void BevView::renderScene(const Beam& beam, const std::vector<const Aperture*>& apertures,
                          const std::vector<double>& boundaries, double extentX, double extentZ) {
    ImGui::TextDisabled("MLC: %.0f mm from source | %.0f mm before isocenter",
                       beam.getSCD(), beam.getSAD() - beam.getSCD());
    const auto available = ImGui::GetContentRegionAvail();
    if (available.x < 1.0f || available.y < 1.0f) return;
    const auto origin = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("##bev_scene", available);
    handleZoom(m_sceneZoom);
    const auto dpi = ImGui::GetIO().DisplayFramebufferScale;
    const auto texture = m_sceneRenderer.render(m_state.patientData.get(), beam, apertures, boundaries,
        extentX, extentZ, std::max(1, static_cast<int>(available.x * dpi.x)),
        std::max(1, static_cast<int>(available.y * dpi.y)), m_sceneZoom);
    if (texture) {
        ImGui::GetWindowDrawList()->AddImage(reinterpret_cast<void*>(static_cast<intptr_t>(texture)),
            origin, ImVec2(origin.x + available.x, origin.y + available.y), ImVec2(0, 1), ImVec2(1, 0));
    } else {
        ImGui::SetCursorScreenPos(origin);
        ImGui::TextUnformatted("Cannot render BEV: invalid source/collimator geometry (see log).");
    }
}

} // namespace optirad
