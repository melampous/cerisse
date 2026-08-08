# Real A/B threshold freeze

Frozen at `2026-07-22T23:49:10Z`, while the Y9000X solver pair was still
running and before either final plotfile existed.

- `analyze_real_ab.py` SHA-256:
  `cfe0d20b4b3cfed54773344aeddab6caf258d71e80c15b27ec33cb29b669d951`
- `verify_real_ab.py` SHA-256:
  `e2bad22578f5580e62b1be2aa0f25b54c0fca97483792d61262890037ac61e0c`

The verifier is allowed to report failure. Its thresholds must not be relaxed
after examining the `real_ab_250_r2` fields.
