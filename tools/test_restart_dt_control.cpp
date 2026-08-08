#include "../src/RestartDtControl.h"

#include <cassert>
#include <cmath>

namespace {

bool close(double lhs, double rhs)
{
  return std::abs(lhs - rhs) <= 1.0e-14 * (1.0 + std::abs(rhs));
}

} // namespace

int main()
{
  using cerisse::time_step::RestartFirstDtFromCflControl;
  using cerisse::time_step::limit_cfl_dt;

  RestartFirstDtFromCflControl control;

  // Fresh/default-off execution never grants an exception.
  assert(!control.pending());
  assert(!control.consume(false));
  control.arm(false);
  assert(!control.pending());

  // A post-regrid estimate neither grows dt nor consumes the restart grant.
  control.arm(true);
  assert(control.pending());
  assert(!control.consume(true));
  assert(control.pending());
  assert(close(limit_cfl_dt(1.0, 0.1, 1.1, true, true), 0.1));

  // The first normal estimate consumes the grant; all later estimates use
  // the legacy 1.1 growth cap again.
  assert(control.consume(false));
  assert(!control.pending());
  assert(!control.consume(false));
  assert(close(limit_cfl_dt(1.0, 0.1, 1.1, false, true), 1.0));
  assert(close(limit_cfl_dt(1.0, 0.1, 1.1, false, false), 0.11));

  // No path increases a CFL estimate that is already below its applicable
  // cap.
  assert(close(limit_cfl_dt(0.05, 0.1, 1.1, false, false), 0.05));
  assert(close(limit_cfl_dt(0.05, 0.1, 1.1, true, false), 0.05));

  return 0;
}
