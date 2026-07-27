FIN P241 TRANSFER TEMPLATE

These files define structure only. They contain no empirical events and
cannot pass the P241 gate.

Production sequence:
1. Provider records raw events and calibration/control/environment files.
2. Registrar replaces template names, computes every SHA-256, and completes
   bundle_manifest.json.
3. Registrar signs bundle_manifest.json with a detached GPG signature.
4. Analyst receives neither the sealed holdout nor registrar private key.
5. Run:
   python3 fin_lab_p241_validator.py BUNDLE      --signature bundle_manifest.json.asc      --trusted-keyring registrar-trustedkeys.gpg      --output admission_certificate.json
