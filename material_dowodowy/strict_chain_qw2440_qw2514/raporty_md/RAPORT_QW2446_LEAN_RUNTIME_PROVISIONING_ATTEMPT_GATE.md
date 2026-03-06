# RAPORT QW-2446: LEAN RUNTIME PROVISIONING ATTEMPT GATE

- Date UTC: 2026-03-05T13:35:59.549400+00:00
- Verdict: **LEAN_RUNTIME_PROVISIONING_ATTEMPT_GATE_PASS_SKIPPED_RUNTIME_ALREADY_AVAILABLE**
- pass_count: `3/9`
- runtime_available_before_attempt: `True`
- provisioning_attempt_executed: `False`
- provisioning_skipped_due_to_existing_runtime: `True`
- curl_download_succeeded: `False`
- dns_resolution_failed: `False`
- runtime_installed_locally: `False`

## Interpretacja
- Runtime Lean był już dostępny przed uruchomieniem tej bramki.
- `runtime_installed_locally=false` oznacza brak nowej instalacji w tej bramce, nie brak runtime.
