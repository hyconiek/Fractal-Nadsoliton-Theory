# Minimal replay

From `fin_physical_bridge_campaign_01/`:

```bash
python -m pytest -q PHYS-002/test_safe_seed.py
python PHYS-002/safe_seed.py --N 2 --g 5.145228719489142
python run_smallN_campaign.py
python PHYS-003/enumerator.py
```

Expected critical checks:
- 4 pytest tests pass;
- N=3 heat-bath k=4 `rho = 0.13143978619564556` (floating replay);
- all 12 C12 sectors are present at N=2,g=0;
- A7 rank 7, centered/symmetric/constant diagonal within floating residuals;
- no calculation uses N>3 in the scientific campaign.

`write_reports.py` only serializes the already computed results/reports; rerun it after `run_smallN_campaign.py` if regenerating documentation.
