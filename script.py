import os

files = []
content_dict = {}

for root, dirs, filenames in os.walk('.'):
    for filename in filenames:
        if filename.endswith('.py'):
            files.append(os.path.join(root, filename))

for file in files:
    with open(file, 'r') as file_content:
        content = file_content.read()
        content_dict[file] = content

with open('QW_STUDIES_SUMMARY_NEW.md', 'w') as summary_file:
    for file, content in content_dict.items():
        if not any(content == content_dict[other_file] for other_file in content_dict if other_file != file):
            summary_file.write(f'## {os.path.basename(file)}\n')
            summary_file.write(f'### Zawartość:\n')
            summary_file.write(content + '\n\n')