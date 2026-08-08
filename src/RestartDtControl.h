#ifndef CERISSE_RESTART_DT_CONTROL_H_
#define CERISSE_RESTART_DT_CONTROL_H_

#include <algorithm>

namespace cerisse::time_step {

// One-shot state used by the level-0 time-step controller.  A post-regrid
// estimate is deliberately not allowed to consume the restart exception:
// post-regrid estimates must remain capped by the pre-regrid time step.
class RestartFirstDtFromCflControl
{
public:
  void arm(bool enabled) noexcept { m_pending = enabled; }

  [[nodiscard]] bool consume(bool post_regrid) noexcept
  {
    if (post_regrid || !m_pending) {
      return false;
    }
    m_pending = false;
    return true;
  }

  [[nodiscard]] bool pending() const noexcept { return m_pending; }

private:
  bool m_pending = false;
};

template <typename Real>
[[nodiscard]] constexpr Real limit_cfl_dt(
    Real cfl_estimate, Real old_dt, Real normal_growth_factor,
    bool post_regrid, bool bypass_normal_growth_limit) noexcept
{
  if (post_regrid) {
    return std::min(cfl_estimate, old_dt);
  }
  if (bypass_normal_growth_limit) {
    return cfl_estimate;
  }
  return std::min(cfl_estimate, normal_growth_factor * old_dt);
}

} // namespace cerisse::time_step

#endif
