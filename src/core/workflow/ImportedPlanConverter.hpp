#pragma once

#include "core/ImportedPlan.hpp"
#include "core/ArcSequence.hpp"
#include "core/Aperture.hpp"
#include "core/workflow/PlanBuilder.hpp"
#include <string>
#include <vector>

namespace optirad {

class Plan;
class Machine;

/// Converts a DICOM RTPLAN (`ImportedPlan`) into optirad planning parameters and delivery sequences.
///
/// Leaf positions are used in IEC BEAM LIMITING DEVICE coordinates. Beam geometry applies the imported collimator
/// rotation so the leaf-aligned coordinates are transformed consistently for dose calculation and visualization.
/// Leaf pairs entirely outside the Y jaws are closed and leaf tips are clipped to the X jaws.
class ImportedPlanConverter {
public:
    /// Arcs of the plan's arc beams, in beam order. The control-point spacing is the finest of the arcs; the
    /// optimization spacing is the multiple of it closest to 10 degrees.
    static VmatSettings toVmatSettings(const ImportedPlan& imported);

    /// optirad machine matching the beam mode of the plan's photon beams (6 MV flattened -> TrueBeamSTx_6X,
    /// 6 MV FFF -> TrueBeamSTx_6FFF); empty when the beams disagree or no model matches.
    static std::string suggestMachine(const ImportedPlan& imported);

    /// Fills fractions, delivery mode, arcs or static beam angles, isocenter, prescription and, when one matches
    /// (see `suggestMachine`), the machine. Bixel width and dose-grid resolution are left to the caller.
    static void applyToPlanConfig(const ImportedPlan& imported, PlanConfig& config);

    /// Delivered VMAT sequence of the plan's arc beams. Slice j of a beam lies between DICOM control points j
    /// and j+1: its aperture is the mean of both, its MU the meterset delivered between them (in optirad weight
    /// units, see `Plan::machineMuPerWeight`). `plan` must have been built with `toVmatSettings` arcs; its sector
    /// width (when it is a VMAT plan with as many arcs) sets the control points' sector indices.
    static ArcSequenceResult toArcSequence(const ImportedPlan& imported, const Plan& plan);

    /// Why `plan` cannot deliver the imported arcs as they are (empty when it can): not VMAT, different arcs
    /// (angles, direction, couch, collimator), isocenter or number of fractions.
    static std::vector<std::string> arcPlanMismatches(const ImportedPlan& imported, const Plan& plan);

    /// Step-and-shoot segments of the plan's static beams, one result per beam in beam order. Each control-point
    /// pair delivering MU becomes a segment with the aperture of its first control point.
    static std::vector<LeafSequenceResult> toLeafSequences(const ImportedPlan& imported, const Plan& plan);

    /// Differences between the plan and what optirad can reproduce with `machine`.
    static std::vector<std::string> validate(const ImportedPlan& imported, const Machine& machine);

    /// Aperture of a control point with leaf pairs outside the Y jaws closed and tips clipped to the X jaws.
    static Aperture collimatedAperture(const ImportedControlPoint& cp, const std::vector<double>& leafBoundaries);
};

} // namespace optirad
