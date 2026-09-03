#pragma once

// AstraDock M20 — Fault detection, isolation & recovery (FDIR).
//
// Physical problem:
//   Sensors lie and actuators degrade: GNSS drops out, star trackers blind,
//   gyros jump bias, thrusters under-deliver, wheels saturate. Classical FDIR
//   detects these from signals the loop ALREADY computes — EKF innovations
//   (NIS chi-square gating), measurement rate-of-change, actuator
//   desired-vs-achieved residuals — isolates the suspect channel, and recovers
//   by exclusion (drop the channel), degradation (coast on predict), or safe
//   mode (hold + null rates). No ML here (M23); every threshold is a chi-square
//   or physics number with a documented basis, never tuned to a plot.
//
// Architecture (never collapse the roles):
//   INJECT (harness-side, truth-adjacent): FaultInjector corrupts measurements
//     / commands per a deterministic schedule. It lives OUTSIDE the flight
//     path — the demo applies it between sensor and filter, exactly where a
//     real fault lives (in the hardware, before the flight software sees it).
//   DETECT (flight-side, estimates + measurements only): InnovationMonitor
//     gates NIS against chi-square 95% bounds with M-of-N persistence (a
//     single NIS spike at 5% false-alarm rate is EXPECTED ~1/20 samples; the
//     M-of-N vote is what makes detection reliable).
//   ISOLATE: per-channel monitors vote; the channel with the sustained worst
//     margin is the suspect. Ambiguity is reported, never forced.
//   RECOVER: exclusion mask drops the channel from updates; safe mode freezes
//     guidance to hold + damps rates; metrics record latency/false alarms.
//
// Units/frames: inherit the monitored signals (NIS dimensionless, residuals in
// signal units). Deterministic: schedule + seed fully determine behavior.

#include "math/vector3.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace astradock::fdir {

namespace detail {

inline void require_positive_fdir(double value, const char* message) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::domain_error(message);
    }
}

}  // namespace detail

// Fault kinds the injector can produce (deterministic schedule entries).
enum class FaultKind {
    none,
    gnss_dropout,       // valid = false over the window
    gnss_bias_jump,     // additive ECI position/velocity bias
    gyro_bias_jump,     // additive body-rate bias (feeds predict)
    star_tracker_drop,  // valid = false over the window
    range_stuck,        // frozen last value (stuck output)
    range_noise_burst,  // inflated noise std over the window
    thruster_degraded,  // achieved thrust scaled by authority factor
    wheel_stuck,        // wheel torque forced to zero (stuck actuator)
};

[[nodiscard]] inline const char* fault_kind_name(FaultKind kind) noexcept {
    switch (kind) {
    case FaultKind::none:
        return "none";
    case FaultKind::gnss_dropout:
        return "gnss_dropout";
    case FaultKind::gyro_bias_jump:
        return "gyro_bias_jump";
    case FaultKind::star_tracker_drop:
        return "star_tracker_drop";
    case FaultKind::range_stuck:
        return "range_stuck";
    case FaultKind::range_noise_burst:
        return "range_noise_burst";
    case FaultKind::thruster_degraded:
        return "thruster_degraded";
    case FaultKind::wheel_stuck:
        return "wheel_stuck";
    case FaultKind::gnss_bias_jump:
        return "gnss_bias_jump";
    }
    return "unknown";
}

// One scheduled fault: kind + window + magnitude vector (meaning per kind:
// bias jump = additive offset; noise burst = extra std; degraded = authority
// factor in [0,1]; dropout/stuck ignore magnitude).
struct FaultSchedule {
    FaultKind kind{FaultKind::none};
    double start_time_s{0.0};
    double end_time_s{0.0};
    math::Vector3 magnitude{};
    double scalar{0.0};

    [[nodiscard]] bool is_active(double time_s) const noexcept {
        return time_s >= start_time_s && time_s <= end_time_s;
    }
};

// Harness-side fault injector: applies the active schedule entries to a
// measurement/command sample. Deterministic; no RNG inside (noise bursts add a
// deterministic bias-shaped offset? NO — bursts scale the CALLER's noise via
// the returned inflated std; see inflated_noise_std()).
struct FaultInjector {
    std::vector<FaultSchedule> schedule{};

    [[nodiscard]] bool dropout_active(FaultKind kind, double time_s) const noexcept {
        for (const auto& entry : schedule) {
            if (entry.kind == kind && entry.is_active(time_s)) {
                return true;
            }
        }
        return false;
    }

    // Additive bias active now for jump kinds (zero vector otherwise).
    [[nodiscard]] math::Vector3 active_bias(FaultKind kind, double time_s) const noexcept {
        for (const auto& entry : schedule) {
            if (entry.kind == kind && entry.is_active(time_s)) {
                return entry.magnitude;
            }
        }
        return math::Vector3{};
    }

    // Authority factor for degradation kinds (1.0 = healthy).
    [[nodiscard]] double authority_factor(FaultKind kind, double time_s) const noexcept {
        for (const auto& entry : schedule) {
            if (entry.kind == kind && entry.is_active(time_s)) {
                return entry.scalar;
            }
        }
        return 1.0;
    }

    // Stuck-output latch: the harness holds the last good value while active.
    [[nodiscard]] bool stuck_active(double time_s) const noexcept {
        return dropout_active(FaultKind::range_stuck, time_s);
    }
};

// M-of-N persistence voter: raises `triggered` only when at least m of the
// last n samples exceed the gate. This is what separates detection (reliable)
// from gating (5% false alarms by construction at the 95% bound).
struct PersistenceVoter {
    std::size_t window_n{5};
    std::size_t threshold_m{3};
    std::vector<int> history{};

    void reset() noexcept {
        history.clear();
    }

    // Returns the current trigger state after pushing `exceed`.
    bool push(bool exceed) {
        if (window_n == 0 || threshold_m == 0 || threshold_m > window_n) {
            throw std::domain_error("Persistence voter needs 0 < m <= n");
        }
        history.push_back(exceed ? 1 : 0);
        if (history.size() > window_n) {
            history.erase(history.begin());
        }
        std::size_t count = 0;
        for (int flag : history) {
            count += static_cast<std::size_t>(flag);
        }
        return history.size() >= threshold_m && count >= threshold_m;
    }
};

// Per-channel NIS innovation monitor: chi-square gate + persistence + margin.
// Threshold basis: 95% upper quantile for the channel DOF (M13 diagnostics:
// 6-DOF GNSS 14.449, 3-DOF star tracker 9.348, 1-DOF range 5.024).
struct InnovationMonitor {
    double gate_threshold{14.4494};
    PersistenceVoter voter{};
    bool triggered{false};
    double worst_margin{0.0};  // max (nis - gate) / gate observed
    double first_trigger_time_s{-1.0};

    void reset() noexcept {
        triggered = false;
        worst_margin = 0.0;
        first_trigger_time_s = -1.0;
        voter.reset();
    }

    // Returns trigger state after processing one NIS sample.
    bool update(double nis, double time_s) {
        if (!std::isfinite(nis) || nis < 0.0 || !std::isfinite(time_s)) {
            throw std::domain_error("NIS monitor inputs must be finite, NIS >= 0");
        }
        const double margin = (nis - gate_threshold) / gate_threshold;
        worst_margin = std::max(worst_margin, margin);
        const bool vote = voter.push(nis > gate_threshold);
        if (vote && !triggered) {
            triggered = true;
            first_trigger_time_s = time_s;
        }
        return triggered;
    }
};

// Actuator tracking residual monitor: |desired - achieved| per axis vs a
// physics bound (e.g. wheel authority after saturation, thruster scale).
// Same M-of-N discipline as the NIS path.
struct ResidualMonitor {
    double bound{0.0};
    PersistenceVoter voter{};
    bool triggered{false};
    double first_trigger_time_s{-1.0};

    void reset() noexcept {
        triggered = false;
        first_trigger_time_s = -1.0;
        voter.reset();
    }

    bool update(double residual, double time_s) {
        if (!std::isfinite(residual) || !std::isfinite(time_s)) {
            throw std::domain_error("Residual monitor inputs must be finite");
        }
        const bool vote = voter.push(residual > bound);
        if (vote && !triggered) {
            triggered = true;
            first_trigger_time_s = time_s;
        }
        return triggered;
    }
};

// Isolation verdict across channels: the triggered channel with the largest
// worst_margin is the prime suspect; multiple triggers => ambiguous flag.
// Never forces certainty: with zero triggers the suspect is "none".
struct IsolationVerdict {
    std::string suspect{"none"};
    bool ambiguous{false};
};

[[nodiscard]] inline IsolationVerdict isolate_fault(
    const std::vector<std::pair<std::string, const InnovationMonitor*>>& channels) {
    IsolationVerdict verdict;
    int triggered_count = 0;
    double best_margin = 0.0;
    for (const auto& [name, monitor] : channels) {
        if (monitor->triggered) {
            ++triggered_count;
            if (monitor->worst_margin > best_margin) {
                best_margin = monitor->worst_margin;
                verdict.suspect = name;
            }
        }
    }
    verdict.ambiguous = triggered_count > 1;
    return verdict;
}

// Recovery actions.
enum class RecoveryAction {
    none,
    exclude_channel,  // drop the suspect measurement channel from updates
    coast_predict,    // dynamics-only predict (no aiding) for a bounded time
    safe_mode,        // freeze guidance to hold + null rates
};

[[nodiscard]] inline const char* recovery_action_name(RecoveryAction action) noexcept {
    switch (action) {
    case RecoveryAction::none:
        return "none";
    case RecoveryAction::exclude_channel:
        return "exclude_channel";
    case RecoveryAction::coast_predict:
        return "coast_predict";
    case RecoveryAction::safe_mode:
        return "safe_mode";
    }
    return "unknown";
}

// Recovery policy: exclusion for single-channel measurement faults, coast for
// ambiguity, safe mode for actuator faults or crush-risk geometry. Pure
// function of verdict + fault class — auditable, no hidden state.
[[nodiscard]] inline RecoveryAction select_recovery(
    const IsolationVerdict& verdict,
    bool actuator_fault,
    bool geometry_critical) {
    if (geometry_critical || actuator_fault) {
        return RecoveryAction::safe_mode;
    }
    if (verdict.suspect == "none") {
        return RecoveryAction::none;
    }
    if (verdict.ambiguous) {
        return RecoveryAction::coast_predict;
    }
    return RecoveryAction::exclude_channel;
}

// Detection performance summary (harness-scored): latency from fault onset to
// first trigger, false-alarm count pre-onset, miss flag.
struct DetectionMetrics {
    double latency_s{-1.0};
    int false_alarms{0};
    bool missed{true};
    bool detected{false};
};

}  // namespace astradock::fdir
