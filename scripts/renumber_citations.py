from pathlib import Path
import re

path = Path(__file__).resolve().parent / '08_build_manuscript.py'
text = path.read_text(encoding='utf-8')

# Keep only references that are actually cited in the manuscript and its
# in-document figure/table legends.  The old numbering was inherited from an
# earlier draft and contained fourteen uncited entries.
cited_old = [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
             19, 20, 24, 25, 31, 33, 34, 37, 38, 39, 40, 41, 42, 43, 44,
             45, 46, 47, 48, 49, 51, 52, 53]
mapping = {old: new for new, old in enumerate(cited_old, start=1)}

def replace_citation(match: re.Match) -> str:
    body = match.group(1)
    parts = [x.strip() for x in body.split(',')]
    raw_numbers = []
    for part in parts:
        if '-' in part:
            a, b = [int(x.strip()) for x in part.split('-', 1)]
            raw_numbers.extend([a, b])
        else:
            raw_numbers.append(int(part))
    if any(x not in mapping for x in raw_numbers):
        return match.group(0)
    out = []
    for part in parts:
        if '-' in part:
            a, b = [int(x.strip()) for x in part.split('-', 1)]
            out.append(f'{mapping[a]}-{mapping[b]}')
        else:
            out.append(str(mapping[int(part)]))
    return '[' + ','.join(out) + ']'

# Only rewrite numeric citation brackets; ordinary Python indexing is untouched.
text = re.sub(r'\[((?:\d+)(?:\s*[-,]\s*\d+)*)\]', replace_citation, text)

old_line = 'REFERENCES = {**ORIGINAL_REFS, **NEW_REFS}'
new_block = '''ALL_REFERENCES = {**ORIGINAL_REFS, **NEW_REFS}
CITED_REFERENCE_NUMBERS = %s
REFERENCE_RENUMBER = {old: new for new, old in enumerate(CITED_REFERENCE_NUMBERS, start=1)}
REFERENCES = {
    REFERENCE_RENUMBER[old]: text
    for old, text in ALL_REFERENCES.items()
    if old in REFERENCE_RENUMBER
}''' % repr(cited_old)
if old_line not in text:
    raise SystemExit('Expected REFERENCES assignment was not found')
text = text.replace(old_line, new_block)
path.write_text(text, encoding='utf-8')
print('Updated', path)
print('Old-to-new mapping:', mapping)
