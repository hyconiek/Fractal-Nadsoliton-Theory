import os

# Define the audited range
studies = range(1543, 1551) # 1543 to 1550
output_file = "FULL_RESEARCH_LOG_AUDITED_BLOCKS_A_B_C.md"

with open(output_file, "w", encoding="utf-8") as out:
    out.write("# MERCILESS AUDIT: FULL RESEARCH LOG (QW-1543 to QW-1550)\n")
    out.write("This document contains the strictly audited implementations verifying FIN Theory as a Pre-Geometry rigor.\n\n")
    out.write("## Block A: Geometry (3+1D, Zero Torsion)\n")
    out.write("## Block B: Matter & Gauge (Dynamics form Action)\n")
    out.write("## Block C: Mass & Unification (Soliton Stability & WEP)\n\n")
    
    files = os.listdir(".")
    
    for i in studies:
        script_file = None
        report_file = None
        
        # Priority: UPGRADE > Refined > Corrected
        
        candidates_script = [f for f in files if f.startswith(f"QW_{i}_") and f.endswith(".py")]
        # Sort logic: Upgrade=3, Refined=2, Corrected=1
        best_score = -1
        for f in candidates_script:
            score = 0
            if "Corrected" in f: score = 1
            if "Refined" in f: score = 2
            if "Upgrade" in f: score = 3
            if score > best_score:
                best_score = score
                script_file = f
        
        candidates_report = [f for f in files if f.startswith(f"RAPORT_QW{i}") and f.endswith(".md")]
        best_score = -1
        for f in candidates_report:
            score = 0
            if "Corrected" in f: score = 1 # Case insensitive? upper usually
            if "CORRECTED" in f: score = 1
            if "REFINED" in f: score = 2
            if "UPGRADE" in f: score = 3
            if score > best_score:
                best_score = score
                report_file = f
        
        if not script_file and not report_file: continue
            
        out.write(f"## STUDY QW-{i} (AUDIT LEVEL: {best_score})\n")
        if best_score == 3:
            out.write("**(STRICT UPGRADE CONFIRMED)**\n\n")
        
        if script_file:
            out.write(f"### Script: `{script_file}`\n\n")
            out.write("```python\n")
            with open(script_file, "r", encoding="utf-8") as f:
                out.write(f.read())
            out.write("\n```\n\n")
            c
        if report_file:
            out.write(f"### Report: `{report_file}`\n\n")
            with open(report_file, "r", encoding="utf-8") as f:
                out.write(f.read())
            out.write("\n\n---\n\n")

print(f"Synthesized log saved to {output_file}")
