# REPLAY
`python -m pytest -q test_safe_seed.py`
`python safe_seed.py --N 2 --g 5.145228719489142`
`python ../run_smallN_campaign.py`

Moduł nie uruchamia kosztownych obliczeń przy imporcie; CLI ogranicza N do 1..4.
