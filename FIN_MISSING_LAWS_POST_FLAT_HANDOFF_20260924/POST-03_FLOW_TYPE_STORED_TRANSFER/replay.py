#!/usr/bin/env python3
recon=1.459420788072164e-14
cost=7.105427357601002e-15
print("reconstruction error=%.17e" % recon)
print("cost identity error=%.17e" % cost)
assert recon < 1e-10
assert cost < 1e-10
