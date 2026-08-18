#pragma once

#include "math/quaternion.hpp"

namespace astradock::attitude {

// Spacecraft attitude representation.
// Represents the instantaneous spatial orientation of the spacecraft body frame
// relative to a reference frame (e.g. ECI or LVLH) as a unit quaternion.
//
// In accordance with M08 scope, this structure represents pure rotational orientation
// and deliberately does not include angular rates, torques, or inertia (scheduled for M09).
struct AttitudeState {
    math::Quaternion orientation{math::Quaternion::identity()};
};

}  // namespace astradock::attitude
