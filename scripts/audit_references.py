from pathlib import Path
from docx import Document
import re

doc = Document(Path(__file__).resolve().parents[1] / 'docs' / 'Manuscript_Inflammation_Research.docx')
paragraphs = [p.text for p in doc.paragraphs]
ref_start = next(i for i, p in enumerate(paragraphs) if p.strip() == 'References')
main_text = '\n'.join(paragraphs[:ref_start])
tokens = []
for block in re.findall(r'\[([^\]]+)\]', main_text):
    for item in re.split(r'[,;]', block):
        item = item.strip()
        if not item:
            continue
        if '-' in item:
            a, b = [int(x) for x in item.split('-', 1)]
            tokens.extend(range(a, b + 1))
        elif item.isdigit():
            tokens.append(int(item))
cited = sorted(set(tokens))
refs = {}
for p in paragraphs[ref_start + 1:]:
    m = re.match(r'^(\d+)\.\s', p)
    if m:
        refs[int(m.group(1))] = p
print('CITED', cited)
print('MISSING', sorted(set(cited) - set(refs)))
print('UNCITED', sorted(set(refs) - set(cited)))
for n in sorted(refs):
    print(f'{n}\t{"CITED" if n in cited else "UNCITED"}\t{refs[n]}')
