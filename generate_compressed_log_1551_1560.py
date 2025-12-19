import os
import re

OUTPUT_FILE = "FULL_LOG_COMPRESSED_PHASE2.md"

STUDIES = [
    ("QW-1551", "QW_1551_Renormalization_Flow.py", "RAPORT_QW1551_RG_FLOW.md"),
    ("QW-1552", "QW_1552_Emergent_Friedmann.py", "RAPORT_QW1552_FRIEDMANN.md"),
    ("QW-1553", "QW_1553_Dark_Energy.py", "RAPORT_QW1553_DARK_ENERGY.md"),
    ("QW-1554", "QW_1554_Dark_Matter_Candidate.py", "RAPORT_QW1554_DARK_MATTER.md"),
    ("QW-1555", "QW_1555_Matter_Geometry_Duality.py", "RAPORT_QW1555_DUALITY.md"),
    ("QW-1556", "QW_1556_Information_Conservation.py", "RAPORT_QW1556_INFORMATION.md"),
    ("QW-1557", "QW_1557_Black_Hole_Information.py", "RAPORT_QW1557_BLACK_HOLE.md"),
    ("QW-1558", "QW_1558_Quantum_Measurement.py", "RAPORT_QW1558_MEASUREMENT.md"),
    ("QW-1560", "QW_1560_TOE_Closure_Test.py", "RAPORT_QW1560_TOE_CLOSURE.md")
]

def clean_code(filename):
    if not os.path.exists(filename):
        return f"[MISSING CODE: {filename}]"
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    clean_lines = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            if not any(keyword in stripped for keyword in ["Upgrade", "QW-", "MERCILESS"]):
                continue
        # Inline comments removal
        if "  # " in line:
            line = line.split("  # ")[0]
        clean_lines.append(line.rstrip())
        
    return "\n".join(clean_lines)

def extract_results(filename):
    if not os.path.exists(filename):
        return f"[MISSING REPORT: {filename}]"
    
    with open(filename, 'r') as f:
        content = f.read()
    
    if "## Results" in content:
        results = content.split("## Results")[1]
    else:
        results = content
        
    return results.strip()

with open(OUTPUT_FILE, 'w') as out:
    out.write("# FULL AUDIT LOG COMPRESSED PHASE 2 (QW-1551 - QW-1560)\n")
    out.write("**Synthesis of Logic & Verification results.**\n\n")
    
    for label, code_file, report_file in STUDIES:
        out.write(f"## {label}\n")
        
        # Code
        code = clean_code(code_file)
        out.write(f"### S:{code_file}\n")
        out.write("```python\n")
        out.write(code)
        out.write("\n```\n\n")
        
        # Results
        res = extract_results(report_file)
        out.write(f"### R:{report_file}\n")
        out.write(res)
        out.write("\n\n" + "="*50 + "\n\n")

print(f"Generated {OUTPUT_FILE}")
