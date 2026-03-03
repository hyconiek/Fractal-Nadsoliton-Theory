import os
import glob
import re

output_file = "phase47_66_code_report.md"

def strip_code(code):
    lines = code.split("\n")
    cleaned = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("import ") or stripped.startswith("from ") or stripped.startswith("logging."):
            continue
        if stripped.startswith("log.") or stripped.startswith("print("):
            continue
        if stripped.startswith("#"):
            continue
        if stripped == "" and (len(cleaned) == 0 or cleaned[-1] == ""):
            continue
        cleaned.append(line)
    return "\n".join(cleaned).strip()

with open(output_file, "w") as out:
    out.write("# Metodologiczny Raport Kodu: Fazy 47-66\n\n")
    out.write("Poniżej znajduje się zagregowany kod ze wszystkich faz (pozbawiony importów, printów, logów i komentarzy w celu ułatwienia oceny matematycznej) oraz wygenerowane wyniki (JSON), jeżeli skrypt je zapisywał.\n\n")

    for i in range(47, 67):
        py_files = glob.glob(f"phase{i}_*.py")
        if not py_files:
            continue
            
        for py_file in py_files:
            with open(py_file, "r") as f:
                code = f.read()
                
            out.write(f"## Faza {i}: `{py_file}`\n")
            out.write("### Kod:\n```python\n")
            out.write(strip_code(code))
            out.write("\n```\n")
            
            # Find associated json
            json_files = glob.glob(f"*v{i}*.json") + glob.glob(f"*{i}*.json")
            json_files = [j for j in json_files if j != output_file] # in case
            
            # Remove duplicates
            json_files = list(set(json_files))
            
            if json_files:
                out.write("\n### Wynik / JSON:\n")
                for jf in json_files:
                    try:
                        with open(jf, "r") as j:
                            out.write("```json\n" + j.read() + "\n```\n")
                    except:
                        pass
            out.write("\n---\n\n")

print(f"Created {output_file}")
