# Deduplication report

- Logical file references across included campaigns: **610**
- Unique payload objects stored: **575**
- Logical bytes before exact-content deduplication: **1906087**
- Unique payload bytes in content store: **1791226**
- Exact duplicate bytes eliminated: **114861**
- Reduction from exact duplicates: **6.03%**

No `.zip` payload is stored. Exactly one `.pyc` is retained because the latest closure's original MANIFEST explicitly requires that hash for exact verifier compatibility; all other cache bytecode is omitted. Earlier package archives are deliberately excluded; their logical content is represented through the campaign maps.

## Campaigns

- `00_phy_full`: 122 logical files; 122 new objects; 0 duplicate references.
- `01_source_geometry`: 68 logical files; 63 new objects; 5 duplicate references.
- `02_dimension_principle_dynamics`: 48 logical files; 43 new objects; 5 duplicate references.
- `03_hdim_source`: 63 logical files; 58 new objects; 5 duplicate references.
- `04_resonance_carrier`: 63 logical files; 60 new objects; 3 duplicate references.
- `05_carrier_source`: 46 logical files; 43 new objects; 3 duplicate references.
- `06_operator_refinement_intermediate`: 21 logical files; 18 new objects; 3 duplicate references.
- `07_operational_carrier`: 25 logical files; 23 new objects; 2 duplicate references.
- `08_ocb_memory_cal_channel`: 25 logical files; 23 new objects; 2 duplicate references.
- `09_hankel_memory`: 24 logical files; 22 new objects; 2 duplicate references.
- `10_hankel2_mem3_cal3_channel3`: 25 logical files; 23 new objects; 2 duplicate references.
- `11_hankel3_mem4_cal4_channel4`: 26 logical files; 25 new objects; 1 duplicate references.
- `12_hankel3b_mem5_cal5_channel5`: 27 logical files; 26 new objects; 1 duplicate references.
- `13_ocb_current_queue_closure`: 27 logical files; 26 new objects; 1 duplicate references.
