from zipfile import ZipFile
from pathlib import Path
from lxml import etree
p=Path(r'docs/Manuscript_Inflammation_Research_final_PMID_manual_EndNote.docx')
with ZipFile(p) as z:
    print('media', [n for n in z.namelist() if n.startswith('word/media/')])
    xml=z.read('word/document.xml')
root=etree.fromstring(xml); ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
text='\n'.join(''.join(t.text or '' for t in para.xpath('.//w:t',namespaces=ns)) for para in root.xpath('.//w:p',namespaces=ns))
print('PMID labels', text.count('PMID:'))
print('\n'.join(text.splitlines()[-20:]))
