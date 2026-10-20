# Bounded output repair retained

Date: 2026-10-20. The first record run reached the log-flow controls but
failed while converting the exact summed log enclosure to text: its
rational denominator exceeded Python's 4300-digit integer string limit.
This was an output-serialization failure, not a failed mathematical
inequality or a passed full-package verification.

Scientific arithmetic still uses the full exact Fraction endpoints.
Only the saved log-form endpoints are now rounded outward to the exact
1/10^30 grid, with their enclosure checked before serialization. The
integer-string safety limit was not disabled and no scientific constant,
inequality, model, input hash or predecessor file was changed. The final
read-only replay must pass after the local output repair.

This is separate from the explicitly recorded analytical refinement
a0 -> a0+a1/N in the proof. It adds no new proof route or physical premise.
