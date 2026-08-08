
Here are the scripts, automatic formatting helpers, yt examples, and validation
tools.

## IBM Mach-4 Circle Gate

`check_ibm_circle_gate.py` enforces the frozen, steady `N=320` GPIBM bow-shock
software regression.  A nonzero exit status is a hard stop for unrelated IBM
development:

```bash
python3 tools/check_ibm_circle_gate.py <circle-history-directory> \
  --output circle_gate.json
```

The default reference is the 2026-07-17 BI-centered constrained shared-GP
LLF-WENO-Z5 result.  This gate does not replace a matched body-fitted physical
validation solution.
