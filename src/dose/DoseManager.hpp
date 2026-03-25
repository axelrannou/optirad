#pragma once

#include "dose/DoseMatrix.hpp"
#include "geometry/Grid.hpp"
#include <memory>
#include <vector>
#include <string>
#include <algorithm>

namespace optirad {

/// A single named dose map with its associated grid.
struct DoseEntry {
    int id = 0;
    std::string name;
    std::shared_ptr<DoseMatrix> dose;
    std::shared_ptr<Grid> grid;
};

/// Manages a collection of dose maps (imported + optimization results).
/// Lives in GuiAppState so all panels can access it.
class DoseManager {
public:
    /// Add a new dose map. Auto-selects it as current. Returns the entry id.
    int addDose(const std::string& name,
                std::shared_ptr<DoseMatrix> dose,
                std::shared_ptr<Grid> grid) {
        DoseEntry entry;
        entry.id = m_nextId++;
        entry.name = name;
        entry.dose = std::move(dose);
        entry.grid = std::move(grid);
        m_entries.push_back(std::move(entry));
        m_selectedIdx = static_cast<int>(m_entries.size()) - 1;
        ++m_version;
        return m_entries.back().id;
    }

    /// Remove dose at given list index.
    void removeDose(int idx) {
        if (idx < 0 || idx >= static_cast<int>(m_entries.size())) return;
        m_entries.erase(m_entries.begin() + idx);
        // Adjust selection
        if (m_entries.empty()) {
            m_selectedIdx = -1;
        } else if (m_selectedIdx >= static_cast<int>(m_entries.size())) {
            m_selectedIdx = static_cast<int>(m_entries.size()) - 1;
        }
        // Adjust comparison index
        if (m_compareIdx == idx) {
            m_compareIdx = -1;
        } else if (m_compareIdx > idx) {
            --m_compareIdx;
        }
        ++m_version;
    }

    /// Get currently selected entry (may be nullptr).
    const DoseEntry* getSelected() const {
        if (m_selectedIdx < 0 || m_selectedIdx >= static_cast<int>(m_entries.size()))
            return nullptr;
        return &m_entries[m_selectedIdx];
    }

    /// Get comparison entry (may be nullptr).
    const DoseEntry* getCompare() const {
        if (m_compareIdx < 0 || m_compareIdx >= static_cast<int>(m_entries.size()))
            return nullptr;
        return &m_entries[m_compareIdx];
    }

    int getSelectedIdx() const { return m_selectedIdx; }
    int getCompareIdx() const { return m_compareIdx; }

    void setSelected(int idx) {
        if (idx >= 0 && idx < static_cast<int>(m_entries.size())) {
            m_selectedIdx = idx;
            ++m_version;
        }
    }

    void setCompare(int idx) {
        if (idx >= -1 && idx < static_cast<int>(m_entries.size())) {
            m_compareIdx = idx;
            ++m_version;
        }
    }

    int count() const { return static_cast<int>(m_entries.size()); }
    int version() const { return m_version; }

    const DoseEntry* getEntry(int idx) const {
        if (idx < 0 || idx >= static_cast<int>(m_entries.size())) return nullptr;
        return &m_entries[idx];
    }

    const std::vector<DoseEntry>& getEntries() const { return m_entries; }

    /// Clear all dose entries.
    void clear() {
        m_entries.clear();
        m_selectedIdx = -1;
        m_compareIdx = -1;
        ++m_version;
    }

    /// Return the next sequential optimization number (for naming).
    int nextOptimizationNumber() const { return m_optimizationCount + 1; }

    /// Increment the optimization counter (call after adding an optimization dose).
    void incrementOptimizationCount() { ++m_optimizationCount; }

private:
    std::vector<DoseEntry> m_entries;
    int m_selectedIdx = -1;
    int m_compareIdx = -1;
    int m_nextId = 1;
    int m_version = 0;
    int m_optimizationCount = 0;
};

} // namespace optirad
