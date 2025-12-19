import os

OUTPUT_FILE = "FULL_LOG_COMPRESSED_PHASE2_UPGRADE.md"

STUDIES = [
    ("QW-1551 (RG Flow)", "QW_1551_Renormalization_Flow.py", "RAPORT_QW1551_RG_FLOW.md"),
    ("QW-1552 (Friedmann)", "QW_1552_Emergent_Friedmann.py", "RAPORT_QW1552_FRIEDMANN.md"),
    ("QW-1553 (Dark Energy)", "QW_1553_Dark_Energy.py", "RAPORT_QW1553_DARK_ENERGY.md"),
    ("QW-1554 (Dark Matter)", "QW_1554_Dark_Matter_Candidate.py", "RAPORT_QW1554_DARK_MATTER.md"),
    ("QW-1548bis (Duality)", "QW_1548bis_Matter_Geometry_Duality.py", "RAPORT_QW1548bis_DUALITY.md"),
    ("QW-1556 (Information)", "QW_1556_Information_Conservation.py", "RAPORT_QW1556_INFORMATION.md"),
    ("QW-1557 (Black Hole)", "QW_1557_Black_Hole_Information.py", "RAPORT_QW1557_BLACK_HOLE.md"),
    ("QW-1558' (Measurement)", "QW_1558_Quantum_Measurement.py", "RAPORT_QW1558_MEASUREMENT_UPGRADE.md"),
    ("QW-1559 (Axioms)", "QW_1559_Minimal_Axioms.py", "RAPORT_QW1559_AXIOMS.md"),
    ("QW-1560 (Classicality)", "QW_1560_Classicality_Emergence.py", "RAPORT_QW1560_CLASSICALITY.md"),
    ("QW-1561 (Closure)", "QW_1561_TOE_Closure_Test.py", "RAPORT_QW1561_TOE_CLOSURE.md")
]

def clean_code(filename):
    if not os.path.exists(filename):
        return f"[MISSING CODE: {filename}]"
    with open(filename, 'r') as f:
        lines = f.readlines()
    clean_lines = []
    for line in lines:
        stripped = line.strip()
        if not stripped: continue
        if stripped.startswith("#"):
            if not any(keyword in stripped for keyword in ["Upgrade", "QW-", "MERCILESS", "Axiom", "Bifurcation"]):
                continue
        clean_lines.append(line.rstrip())
    return "\n".join(clean_lines)

def extract_results(filename):
    if not os.path.exists(filename):
        return f"[MISSING REPORT: {filename}]"
    with open(filename, 'r', encoding="utf-8") as f:
        content = f.read()
    
    # Extract the whole report but maybe highlight results
    return content.strip()

with open(OUTPUT_FILE, 'w', encoding="utf-8") as out:
    out.write("# FULL AUDIT LOG COMPRESSED PHASE 2 UPGRADE (QW-1551 - QW-1561)\n")
    out.write("**Strict Audit Version - Zero Assumption Physics.**\n\n")
    
    for label, code_file, report_file in STUDIES:
        out.write(f"## {label}\n")
        code = clean_code(code_file)
        out.write(f"### S:{code_file}\n")
        out.write("```python\n")
        out.write(code)
        out.write("\n```\n\n")
        res = extract_results(report_file)
        out.write(f"### R:{report_file}\n\n")
        out.write(res)
        out.write("\n\n" + "="*80 + "\n\n")

print(f"Generated {OUTPUT_FILE}")
