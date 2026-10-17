# Reproduction instructions

Run `python verify_dh.py` in the campaign root. The script performs exact integer group/orbit geometry and graph predecessor-distance checks, and exact fraction tests of the complement inequalities. These are bounded controls; they DO NOT prove the uniform Stirling/lattice lower mass or MT E barrier. Run `python check_package.py` after creating the ZIP, or simply examine `MANIFEST.sha256` and archives. Source file hashes in `INPUTS.sha256` are of accessible local bytes, not missing actual current Git snapshots.
