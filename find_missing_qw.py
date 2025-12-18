
import os
import re

def get_files_qw():
    qw_pattern = re.compile(r'QW-V?\d+')
    files = os.listdir('.')
    qw_files = set()
    for f in files:
        matches = qw_pattern.findall(f)
        for m in matches:
            num = int(re.search(r'\d+', m).group())
            if num >= 400:
                qw_files.add(m)
    return sorted(list(qw_files), key=lambda x: (int(re.search(r'\d+', x).group()), 'V' in x))

def get_sum_qw():
    qw_pattern = re.compile(r'QW-V?\d+')
    with open('gemini_sum.md', 'r') as f:
        content = f.read()
    return set(qw_pattern.findall(content))

all_qw = get_files_qw()
sum_qw = get_sum_qw()

missing = [qw for qw in all_qw if qw not in sum_qw]
print(f"Missing QW entries (400+): {missing}")
print(f"Total all: {len(all_qw)}, Total in sum: {len(sum_qw)}, Total missing: {len(missing)}")
