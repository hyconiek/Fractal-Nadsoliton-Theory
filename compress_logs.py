import os
import re
import json

SOURCE_DIR = "/home/krzysiek/Pobrane/TOE/edison"
OUTPUT_FILE = os.path.join(SOURCE_DIR, "FULL_LOG_COMPRESSED_QW420_1200.md")
MIN_QW = 420
MAX_QW = 1200

def get_qw_range(filename):
    """
    Extracts QW numbers from filename.
    Handles:
    - QW-123
    - QW_123
    - QW-100_to_QW-200
    - QW-100-QW-200
    Returns a set of numbers covered by the file.
    """
    # Find all QW numbers
    matches = re.findall(r'QW[-_]?(\d+)', filename, re.IGNORECASE)
    if not matches:
        return set()
    
    numbers = [int(m) for m in matches]
    
    # Check if it looks like a range (e.g. "to", "through" or just distinct numbers implying a batch)
    # If the filename has "to" or "through" between numbers, assume range.
    # Simple heuristic: if we have 2 numbers and they are far apart, maybe it's a range?
    # But usually files are named explicitly "QW-A_to_QW-B".
    
    if len(numbers) == 2 and ("to" in filename.lower() or "through" in filename.lower() or "-" in filename):
        # Assume range
        start, end = min(numbers), max(numbers)
        if end - start < 1000: # Sanity check for massive ranges
            return set(range(start, end + 1))
            
    return set(numbers)

def is_relevant(filename):
    if not (filename.endswith('.py') or filename.endswith('.md') or filename.endswith('.json')):
        return False
    if filename.startswith("FULL_LOG_COMPRESSED"):
        return False
        
    qw_nums = get_qw_range(filename)
    if not qw_nums:
        return False
        
    # Check intersection
    for num in qw_nums:
        if MIN_QW <= num <= MAX_QW:
            return True
    return False

def compress_python(content):
    lines = content.splitlines()
    compressed = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("import ") or stripped.startswith("from "):
            continue
        # Remove heavy banner comments potentially
        if stripped.startswith("# ==") or stripped.startswith("# --"):
            continue
        compressed.append(line)
    return "\n".join(compressed)

def compress_markdown(content):
    lines = content.splitlines()
    compressed = []
    for line in lines:
        if not line.strip():
            continue
        compressed.append(line)
    return "\n".join(compressed)

def compress_json(content):
    try:
        data = json.loads(content)
        # return compact json
        return json.dumps(data, separators=(',', ':'))
    except:
        return content

def main():
    files = [f for f in os.listdir(SOURCE_DIR) if is_relevant(f)]
    # Sort files - try to sort by first QW number found
    def sort_key(f):
        nums = list(get_qw_range(f))
        return nums[0] if nums else 999999
        
    files.sort(key=sort_key)
    
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as outfile:
        outfile.write(f"# COMPRESSED LOG: QW-{MIN_QW} to QW-{MAX_QW}\n\n")
        
        for f in files:
            path = os.path.join(SOURCE_DIR, f)
            try:
                with open(path, 'r', encoding='utf-8', errors='replace') as infile:
                    content = infile.read()
                
                outfile.write(f"## FILE: {f}\n")
                
                if f.endswith('.py'):
                    c_content = compress_python(content)
                    outfile.write("```python\n")
                    outfile.write(c_content)
                    outfile.write("\n```\n\n")
                elif f.endswith('.md'):
                    c_content = compress_markdown(content)
                    outfile.write(c_content)
                    outfile.write("\n\n")
                elif f.endswith('.json'):
                    c_content = compress_json(content)
                    outfile.write("```json\n")
                    outfile.write(c_content)
                    outfile.write("\n```\n\n")
                    
            except Exception as e:
                outfile.write(f"Error reading {f}: {e}\n\n")

    print(f"Created {OUTPUT_FILE} with {len(files)} files processed.")

if __name__ == "__main__":
    main()
