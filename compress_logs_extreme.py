import os
import re
import json

SOURCE_DIR = "/home/krzysiek/Pobrane/TOE/edison"
OUTPUT_FILE = os.path.join(SOURCE_DIR, "FULL_LOG_COMPRESSED_EXTREME_QW420_1200_FINAL_V3.md")
MIN_QW = 420
MAX_QW = 1200

# STRICT LINE LIMITS
MAX_PY_LINES = 50
MAX_MD_LINES = 30
MAX_JSON_LINES = 20

def get_qw_numbers(filename):
    matches = re.findall(r'QW[-_]?(\d+)', filename, re.IGNORECASE)
    if not matches:
        return set()
    nums = [int(m) for m in matches]
    if len(nums) == 2 and ("to" in filename.lower() or "through" in filename.lower() or "-" in filename):
        start, end = min(nums), max(nums)
        if end - start < 1000:
            return set(range(start, end + 1))
    return set(nums)

def is_relevant_file(filename):
    if filename.startswith("FULL_LOG_COMPRESSED"): return False
    if not (filename.endswith('.py') or filename.endswith('.md') or filename.endswith('.json')): return False
    
    nums = get_qw_numbers(filename)
    if not nums: return False
    return any(MIN_QW <= n <= MAX_QW for n in nums)

def extract_methodology_from_py_code(content):
    """
    Extracts purely executable methodology logic with VERY STRICT heuristics.
    """
    lines = content.splitlines()
    extracted = []
    in_docstring_double = False
    in_docstring_single = False

    for line in lines:
        stripped = line.strip()
        if not stripped: continue
        
        # Docstring Logic
        cnt_double = stripped.count('"""')
        if cnt_double % 2 == 1:
            in_docstring_double = not in_docstring_double
            continue 
        
        cnt_single = stripped.count("'''")
        if cnt_single % 2 == 1:
             in_docstring_single = not in_docstring_single
             continue

        if in_docstring_double or in_docstring_single:
            continue
            
        if (stripped.startswith('"""') and stripped.endswith('"""')) or \
           (stripped.startswith("'''") and stripped.endswith("'''")):
            continue

        # Exclusions
        if stripped.startswith("#"): continue
        if stripped.startswith("import ") or stripped.startswith("from "): continue
        if any(x in stripped for x in ["print(", "plt.", "ax.", "figure(", "show()", "savefig(", "title(", "xlabel(", "ylabel("]): continue
        
        # Markdown Text Avoidance
        # Skip lines starting with markdown indicators usually found in my descriptions
        if stripped.startswith("*") or stripped.startswith("- ") or stripped.startswith("|") or stripped.startswith(">"):
             continue

        # STRICT Logic Heuristics
        is_logic = False
        
        # 1. Math/Assignment
        if "=" in stripped:
             # Exclude if line looks like "Word = Word" text, unless it has code symbols
             if re.search(r'[+\-*/%<>&|^]', stripped) or "[" in stripped or "(" in stripped or "np." in stripped or "self." in stripped:
                 is_logic = True
             elif re.match(r'^[a-zA-Z0-9_]+\s*=\s*', stripped): # simple assignment: x = 5
                 is_logic = True
        
        # 2. Control Flow Keywords
        elif re.search(r'\b(if|for|while|with|def|class|return|yield|try|except|elif|else)\b', stripped):
             is_logic = True
             
        # 3. Function Calls (Start of line or after space)
        # e.g. "func(a)" or "obj.method()"
        elif re.match(r'^[a-zA-Z0-9_.]+\s*\(', stripped):
             is_logic = True
             
        # 4. Data Structures
        elif stripped.startswith("[") or stripped.startswith("{"):
             is_logic = True

        if is_logic:
            extracted.append(line)

    if len(extracted) > MAX_PY_LINES:
        return "\n".join(extracted[:MAX_PY_LINES]) + "\n... [TRUNCATED LOGIC]"
    return "\n".join(extracted)

def extract_results_from_md(content):
    lines = content.splitlines()
    extracted = []
    keep_block = False
    
    for line in lines:
        stripped = line.strip()
        if not stripped: continue
        
        if stripped.startswith("#"):
            upper = stripped.upper()
            if any(k in upper for k in ["RESULT", "VERDICT", "CONCLUSION", "STATUS", "SUMMARY", "FINDING", "RAPORT", "ANALIZA", "ODKRYCIE", "QW-"]):
                keep_block = True
                extracted.append(line)
            else:
                keep_block = False
            continue
            
        if stripped.startswith("|"):
            extracted.append(line)
            continue
            
        if "**" in stripped:
            if any(k in stripped.upper() for k in ["PASSED", "FAILED", "SUKCES", "PORAŻKA", "VERIFIED", "CONSISTENT", "DERIVED", "ERROR", "%", "RESULT"]):
                extracted.append(line)
            continue
            
        if keep_block and len(stripped) < 200:
            extracted.append(line)

    if len(extracted) > MAX_MD_LINES:
         return "\n".join(extracted[:MAX_MD_LINES]) + "\n... [TRUNCATED RESULTS]"
    return "\n".join(extracted)

def extract_results_from_json(content):
    try:
        data = json.loads(content)
        valid_keys = ["result", "verdict", "status", "error", "conclusion", "summary", "finding", "qw_id", "test_type", "val", "value", "parameters"]
        lines = []
        def recursive_search(obj, depth=0):
            if depth > 3: return
            if len(lines) > MAX_JSON_LINES: return
            
            if isinstance(obj, dict):
                for k, v in obj.items():
                    if any(vk in k.lower() for vk in valid_keys):
                         lines.append(f"{'  '*depth}- **{k}**: {v}")
                    elif isinstance(v, (dict, list)):
                        recursive_search(v, depth+1)
            elif isinstance(obj, list):
                for item in obj[:3]: 
                    recursive_search(item, depth+1)

        recursive_search(data)
        if len(lines) > MAX_JSON_LINES:
            return "\n".join(lines[:MAX_JSON_LINES]) + "\n... [TRUNCATED JSON]"
        return "\n".join(lines)
    except:
        return "[Invalid JSON]"

def main():
    all_files = [f for f in os.listdir(SOURCE_DIR) if is_relevant_file(f)]
    sorted_files = sorted(all_files, key=lambda x: list(get_qw_numbers(x))[0] if get_qw_numbers(x) else 999999)
    
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as outfile:
        outfile.write(f"# DATA LOG (CODE LOGIC + RESULTS) < 1MB: QW-{MIN_QW} to QW-{MAX_QW}\n")
        outfile.write(f"Extraction Strict Limits: Max {MAX_PY_LINES} lines code, {MAX_MD_LINES} lines results per file.\n\n")
        
        for f in sorted_files:
            path = os.path.join(SOURCE_DIR, f)
            try:
                with open(path, 'r', encoding='utf-8', errors='replace') as infile:
                    content = infile.read()
                
                if f.endswith('.py'):
                    c_content = extract_methodology_from_py_code(content)
                    if c_content.strip():
                        outfile.write(f"## {f} [PY: LOGIC]\n")
                        outfile.write("```python\n")
                        outfile.write(c_content)
                        outfile.write("\n```\n")
                elif f.endswith('.md'):
                    c_content = extract_results_from_md(content)
                    if c_content.strip():
                        outfile.write(f"## {f} [MD: RESULTS]\n")
                        outfile.write(c_content)
                        outfile.write("\n\n---\n\n")
                elif f.endswith('.json'):
                    c_content = extract_results_from_json(content)
                    if c_content.strip():
                        outfile.write(f"## {f} [JSON: DATA]\n")
                        outfile.write(c_content)
                        outfile.write("\n\n---\n\n")
                    
            except Exception as e:
                outfile.write(f"Error reading {f}: {e}\n")

    print(f"Created {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
