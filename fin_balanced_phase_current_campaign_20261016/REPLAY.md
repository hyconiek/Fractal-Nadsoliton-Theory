# Replay

Run `python verify_pc.py` (Python standard library, Fraction exact arithmetic); it writes `verification.json`. The tests check rational weight/current/Hall margins, 36 Fourier directions, a strict negative Hall-control, and one exact iid q12 multinomial balance fixture. They are not a proof for all N nor a replay of unavailable source PC03/PC04.

Then use `python check_package.py <zip_path>` after package creation to validate hashes and ZIP CRC. `MANIFEST.sha256` covers emitted new campaign outputs, excluding itself and the archival ZIP wrapper, without altering original source files.
