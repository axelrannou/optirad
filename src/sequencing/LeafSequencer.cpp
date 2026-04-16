#include "LeafSequencer.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <map>
#include <set>

namespace optirad {

LeafSequenceResult LeafSequencer::sequenceBeam(
    const FluenceMap& fluence,
    const MachineGeometry& mlc,
    const LeafSequencerOptions& opts) {

    LeafSequenceResult result;

    int nRows = fluence.getNumRows();  // Z direction (perpendicular to leaf motion)
    int nCols = fluence.getNumCols();  // X direction (parallel to leaf motion)
    if (nRows == 0 || nCols == 0) return result;

    double calFac = fluence.getMaxFluence();
    if (calFac <= 0.0) return result;

    int numLevels = std::max(1, opts.numLevels);
    result.calFac = calFac;
    result.numLevels = numLevels;

    double bw = fluence.getBixelWidth();
    double originX = fluence.getOriginX();
    double originZ = fluence.getOriginZ();
    double resolution = opts.leafPositionResolution;

    result.originX = originX;

    // ── Step 0: Build full leaf pair boundary array in BEV coordinates ──
    int numPairs = (mlc.numLeaves > 0) ? mlc.numLeaves / 2 : 0;
    std::vector<double> fullLeafBounds; // size = numPairs + 1

    if (numPairs > 0 && !mlc.leafWidths.empty()) {
        fullLeafBounds.push_back(0.0);
        if (mlc.leafWidths.size() == 2) {
            double innerW = mlc.leafWidths[0];
            double outerW = mlc.leafWidths[1];
            int innerPairs = (mlc.numInnerPairs > 0) ? mlc.numInnerPairs : 40;
            int outerPerSide = (numPairs - innerPairs) / 2;
            for (int i = 0; i < outerPerSide; ++i)
                fullLeafBounds.push_back(fullLeafBounds.back() + outerW);
            for (int i = 0; i < innerPairs; ++i)
                fullLeafBounds.push_back(fullLeafBounds.back() + innerW);
            for (int i = 0; i < outerPerSide; ++i)
                fullLeafBounds.push_back(fullLeafBounds.back() + outerW);
        } else {
            double w = mlc.leafWidths[0];
            for (int i = 0; i < numPairs; ++i)
                fullLeafBounds.push_back(fullLeafBounds.back() + w);
        }
        double totalSpan = fullLeafBounds.back();
        double offset = totalSpan / 2.0;
        for (auto& b : fullLeafBounds)
            b -= offset;
    }

    // ── Step 1: Determine field leaf pairs (all that overlap the fluence grid Z extent) ──
    double fluenceZMin = originZ - bw * 0.5;
    double fluenceZMax = originZ + (nRows - 1) * bw + bw * 0.5;

    std::vector<int> fieldLeafPairs;

    if (!fullLeafBounds.empty()) {
        for (int lp = 0; lp < numPairs; ++lp) {
            double lpBot = fullLeafBounds[lp];
            double lpTop = fullLeafBounds[lp + 1];
            // Check overlap with fluence extent
            if (lpTop > fluenceZMin && lpBot < fluenceZMax) {
                fieldLeafPairs.push_back(lp);
            }
        }
    }

    int numFieldLP = static_cast<int>(fieldLeafPairs.size());
    bool hasLeafInfo = (numFieldLP > 0);

    // Fallback: if no MLC info, treat each bixel row as its own "leaf pair"
    if (!hasLeafInfo) {
        numFieldLP = nRows;
        fieldLeafPairs.resize(nRows);
        std::iota(fieldLeafPairs.begin(), fieldLeafPairs.end(), 0);
    }

    // Build leaf pair Z boundaries for field leaf pairs
    if (hasLeafInfo) {
        result.leafPairBoundariesZ.resize(numFieldLP + 1);
        for (int i = 0; i < numFieldLP; ++i)
            result.leafPairBoundariesZ[i] = fullLeafBounds[fieldLeafPairs[i]];
        result.leafPairBoundariesZ[numFieldLP] =
            fullLeafBounds[fieldLeafPairs.back() + 1];
    }

    // Map leaf pair index → local index (0..numFieldLP-1)
    std::map<int, int> lpToLocal;
    for (int i = 0; i < numFieldLP; ++i)
        lpToLocal[fieldLeafPairs[i]] = i;

    // ── Step 2: Resample fluence from bixel grid to leaf-pair grid ──
    // Each leaf pair gets an area-weighted average of overlapping bixel rows.
    std::vector<double> resampledFluence(static_cast<size_t>(numFieldLP) * nCols, 0.0);

    if (hasLeafInfo) {
        for (int li = 0; li < numFieldLP; ++li) {
            double lpBot = result.leafPairBoundariesZ[li];
            double lpTop = result.leafPairBoundariesZ[li + 1];
            double lpWidth = lpTop - lpBot;
            if (lpWidth <= 0.0) continue;

            for (int r = 0; r < nRows; ++r) {
                double rowCenter = originZ + r * bw;
                double rowBot = rowCenter - bw * 0.5;
                double rowTop = rowCenter + bw * 0.5;
                double overlap = std::max(0.0, std::min(lpTop, rowTop) - std::max(lpBot, rowBot));
                if (overlap <= 0.0) continue;

                double weight = overlap / lpWidth;
                for (int c = 0; c < nCols; ++c) {
                    resampledFluence[static_cast<size_t>(li) * nCols + c] +=
                        weight * fluence.getValue(r, c);
                }
            }
        }
    } else {
        // No MLC: 1:1 copy
        for (int r = 0; r < nRows; ++r)
            for (int c = 0; c < nCols; ++c)
                resampledFluence[static_cast<size_t>(r) * nCols + c] = fluence.getValue(r, c);
    }

    // ── Step 3: Quantize resampled fluence on leaf-pair grid ──
    std::vector<int> D_lp(static_cast<size_t>(numFieldLP) * nCols, 0);
    int maxLevel = 0;
    for (int li = 0; li < numFieldLP; ++li) {
        for (int c = 0; c < nCols; ++c) {
            double val = resampledFluence[static_cast<size_t>(li) * nCols + c];
            int level = static_cast<int>(std::round(val / calFac * numLevels));
            level = std::clamp(level, 0, numLevels);
            D_lp[static_cast<size_t>(li) * nCols + c] = level;
            maxLevel = std::max(maxLevel, level);
        }
    }

    // ── Step 4: Back-project quantized fluence to bixel grid for DeliverableDoseCalculator ──
    auto rowToLeafPair = fluence.mapToLeafPairs(mlc);
    std::vector<int> D_bixel(static_cast<size_t>(nRows) * nCols, 0);
    for (int r = 0; r < nRows; ++r) {
        int lp = rowToLeafPair[r];
        auto it = lpToLocal.find(lp);
        if (it == lpToLocal.end()) continue;
        int li = it->second;
        for (int c = 0; c < nCols; ++c) {
            D_bixel[static_cast<size_t>(r) * nCols + c] =
                D_lp[static_cast<size_t>(li) * nCols + c];
        }
    }
    result.quantizedFluence = D_bixel;

    if (maxLevel == 0) return result;

    // ── Step 5: Level extraction on leaf-pair grid ──
    auto snapToResolution = [resolution](double pos) -> double {
        if (resolution <= 0.0) return pos;
        return std::round(pos / resolution) * resolution;
    };

    auto buildAperture = [&](const std::vector<std::vector<bool>>& shape,
                             double muWeight, int segIdx) -> Aperture {
        Aperture aperture;
        aperture.bankA.resize(numFieldLP);
        aperture.bankB.resize(numFieldLP);
        aperture.weight = muWeight;
        aperture.segmentIndex = segIdx;

        for (int li = 0; li < numFieldLP; ++li) {
            int leftOpen = -1, rightOpen = -1;
            for (int c = 0; c < nCols; ++c) {
                if (shape[li][c]) {
                    if (leftOpen < 0) leftOpen = c;
                    rightOpen = c;
                }
            }

            if (leftOpen < 0) {
                // Closed: find center from quantized fluence
                int firstCol = -1, lastCol = -1;
                for (int c = 0; c < nCols; ++c) {
                    if (D_lp[static_cast<size_t>(li) * nCols + c] > 0) {
                        if (firstCol < 0) firstCol = c;
                        lastCol = c;
                    }
                }
                if (firstCol >= 0) {
                    double center = originX + (firstCol + lastCol) * 0.5 * bw;
                    aperture.bankA[li] = snapToResolution(center);
                    aperture.bankB[li] = snapToResolution(center);
                } else {
                    aperture.bankA[li] = 0.0;
                    aperture.bankB[li] = 0.0;
                }
            } else {
                double leftEdge = originX + leftOpen * bw - bw * 0.5;
                double rightEdge = originX + rightOpen * bw + bw * 0.5;
                aperture.bankA[li] = snapToResolution(leftEdge);
                aperture.bankB[li] = snapToResolution(rightEdge);
            }
        }
        return aperture;
    };

    struct SegmentInfo {
        Aperture aperture;
        std::vector<std::vector<bool>> shape;
        int levelWeight;
    };
    std::vector<SegmentInfo> allSegments;

    std::vector<std::vector<bool>> prevShape;
    int currentShapeWeight = 0;

    for (int level = 1; level <= maxLevel; ++level) {
        std::vector<std::vector<bool>> shape(numFieldLP, std::vector<bool>(nCols, false));
        for (int li = 0; li < numFieldLP; ++li)
            for (int c = 0; c < nCols; ++c)
                shape[li][c] = (D_lp[static_cast<size_t>(li) * nCols + c] >= level);

        currentShapeWeight++;
        bool isDifferent = (level == 1) || (shape != prevShape);

        if (isDifferent && level > 1 && currentShapeWeight > 1) {
            int w = currentShapeWeight - 1;
            double wt = static_cast<double>(w) / numLevels * calFac;
            auto ap = buildAperture(prevShape, wt, static_cast<int>(allSegments.size()));
            allSegments.push_back({std::move(ap), prevShape, w});
            currentShapeWeight = 1;
        }

        prevShape = shape;

        if (level == maxLevel) {
            int w = currentShapeWeight;
            double wt = static_cast<double>(w) / numLevels * calFac;
            auto ap = buildAperture(shape, wt, static_cast<int>(allSegments.size()));
            allSegments.push_back({std::move(ap), shape, w});
        }
    }

    // ── Step 6: Post-hoc merge small segments ──
    if (opts.minSegmentMU > 0.0 && allSegments.size() > 1) {
        bool changed = true;
        while (changed) {
            changed = false;
            for (size_t i = 0; i < allSegments.size(); ++i) {
                if (allSegments[i].aperture.weight < opts.minSegmentMU) {
                    if (allSegments.size() == 1) break;

                    if (i + 1 < allSegments.size()) {
                        allSegments[i].aperture.weight += allSegments[i + 1].aperture.weight;
                        allSegments[i].levelWeight += allSegments[i + 1].levelWeight;
                        allSegments.erase(allSegments.begin() + static_cast<long>(i + 1));
                    } else {
                        allSegments[i - 1].aperture.weight += allSegments[i].aperture.weight;
                        allSegments[i - 1].levelWeight += allSegments[i].levelWeight;
                        allSegments.erase(allSegments.begin() + static_cast<long>(i));
                    }
                    changed = true;
                    break;
                }
            }
        }
    }

    // ── Step 7: Build result segments and reconstruct filtered leaf-pair quantized fluence ──
    std::vector<int> filteredD_lp(D_lp.size(), 0);

    for (size_t i = 0; i < allSegments.size(); ++i) {
        auto& seg = allSegments[i];
        seg.aperture.segmentIndex = static_cast<int>(i);
        result.segments.push_back(std::move(seg.aperture));

        for (int li = 0; li < numFieldLP; ++li)
            for (int c = 0; c < nCols; ++c)
                if (seg.shape[li][c])
                    filteredD_lp[static_cast<size_t>(li) * nCols + c] += seg.levelWeight;
    }

    // Re-project filtered leaf-pair fluence back to bixel grid for quantizedFluence
    for (int r = 0; r < nRows; ++r) {
        int lp = rowToLeafPair[r];
        auto it = lpToLocal.find(lp);
        if (it == lpToLocal.end()) continue;
        int li = it->second;
        for (int c = 0; c < nCols; ++c) {
            result.quantizedFluence[static_cast<size_t>(r) * nCols + c] =
                filteredD_lp[static_cast<size_t>(li) * nCols + c];
        }
    }

    // Store leaf-pair-resolution deliverable fluence for BEV rendering
    result.leafPairFluenceCols = nCols;
    result.leafPairFluence.resize(static_cast<size_t>(numFieldLP) * nCols);
    for (int li = 0; li < numFieldLP; ++li)
        for (int c = 0; c < nCols; ++c)
            result.leafPairFluence[static_cast<size_t>(li) * nCols + c] =
                static_cast<double>(filteredD_lp[static_cast<size_t>(li) * nCols + c])
                / numLevels * calFac;

    // ── Step 8: Compute total MU ──
    result.totalMU = 0.0;
    for (const auto& seg : result.segments)
        result.totalMU += seg.weight;

    // ── Step 9: Compute fluence fidelity (Pearson correlation on bixel grid) ──
    const auto& original = fluence.getData();
    size_t n = original.size();
    std::vector<double> deliverableFluence(n);
    for (size_t i = 0; i < n; ++i)
        deliverableFluence[i] = static_cast<double>(result.quantizedFluence[i]) / numLevels * calFac;

    double sumOrig = 0, sumDel = 0;
    for (size_t i = 0; i < n; ++i) { sumOrig += original[i]; sumDel += deliverableFluence[i]; }
    double meanOrig = sumOrig / n;
    double meanDel = sumDel / n;

    double cov = 0, varOrig = 0, varDel = 0;
    for (size_t i = 0; i < n; ++i) {
        double dO = original[i] - meanOrig;
        double dD = deliverableFluence[i] - meanDel;
        cov += dO * dD;
        varOrig += dO * dO;
        varDel += dD * dD;
    }

    if (varOrig > 0 && varDel > 0) {
        result.fluenceFidelity = cov / std::sqrt(varOrig * varDel);
    } else {
        result.fluenceFidelity = (varOrig == 0 && varDel == 0) ? 1.0 : 0.0;
    }

    return result;
}

} // namespace optirad
