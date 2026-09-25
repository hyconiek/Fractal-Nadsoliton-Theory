#!/usr/bin/env python3
"""Independently evaluate the colored-tree entropy limit from its recurrence."""
import mpmath as mp

mp.mp.dps = 80
for colors in (4, 8, 12):
    log_count = mp.log(colors)
    correction = mp.mpf('0')
    for level in range(100):
        small = mp.log1p(mp.exp(-log_count))
        correction += (small - mp.log(2)) / mp.power(2, level + 1)
        log_count = 2 * log_count + small - mp.log(2)
    entropy = mp.log(colors) + correction
    print(colors, mp.nstr(entropy, 30))
