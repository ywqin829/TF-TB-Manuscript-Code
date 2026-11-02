from pathlib import Path
import re, json
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_SECTION
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

root = Path(r'F:/TF_TB_2026')
src = root/'docs/markdown/Manuscript_Inflammation_Research.md'
supp_src = root/'docs/markdown/Supplementary_Information.md'
out = root/'docs/Manuscript_Inflammation_Research_final_PMID_inline_clarified.docx'
text = src.read_text(encoding='utf-8-sig')
lines = text.splitlines()
supp_text = supp_src.read_text(encoding='utf-8-sig')
supp_lines = supp_text.splitlines()

# PMID values deliberately limited to IDs verified in the project audit.  References
# without a PubMed record are labelled explicitly; other unresolved records are
# marked as not yet verified rather than being described as missing references.
pmids = {
  9:'22872737', 10:'29050771', 11:'25703559',
  14:'27734027', 16:'20725040', 17:'29921861', 18:'35841871', 19:'10761931',
  20:'24113868', 37:'42841723', 38:'37690670', 39:'31340985', 40:'36699385', 41:'23323831',
  42:'No PubMed PMID (JSS software paper)', 43:'21414208', 44:'17099194', 45:'31126321', 46:'31510660',
  47:'35672358', 48:'42343006', 49:'42662901', 50:'17921069', 51:'27304615', 52:'25605792',
  53:'23056259', 2:'No PubMed PMID (WHO report)', 3:'29996783', 4:'20808293', 5:'25196020', 6:'27368311', 7:'25765201', 12:'26409285'
}

UNVERIFIED = 'PMID not yet verified'

# Parse references and separate figure legends from body.
ref_start = next(i for i,l in enumerate(lines) if l.strip()=='# References')
refs=[]; i=ref_start+1
while i<len(lines):
    m=re.match(r'^\*\*(\d+)\.\*\*\s*(.*)$', lines[i].strip())
    if m:
        n=int(m.group(1)); body=m.group(2).strip(); i+=1
        cont=[]
        while i<len(lines) and lines[i].strip() and not re.match(r'^\*\*\d+\.\*\*', lines[i].strip()):
            # stop if a new section begins
            if lines[i].startswith('# '): break
            cont.append(lines[i].strip()); i+=1
        if cont: body += ' ' + ' '.join(cont)
        refs.append((n,body))
    else: i+=1

# Figure legends in the main manuscript.
legends=[]
for l in lines:
    s=l.strip()
    if re.match(r'^(\*\*)?Fig\.\s*[0-9S]', s) or s.startswith('Figure S') or s=='Supplementary figure legends':
        s=re.sub(r'^\*\*([^*]+)\*\*', r'\1', s)
        legends.append(s)

# Supplementary figure legends are stored in a separate Markdown file.  Include the
# five actual legends, but omit section headings and figure-link lines.
for l in supp_lines:
    s=l.strip()
    if re.match(r'^Figure S[1-5]\.', s) or re.match(r'^Fig\. S[1-5]\s', s):
        legends.append(s)

# Build a clean, text-only manuscript. Figure link/image lines are intentionally omitted.
def clean_inline(s):
    s=re.sub(r'!\[[^\]]*\]\([^)]*\)', '', s)
    s=re.sub(r'\[([^\]]+)\]\([^)]*\)', r'\1', s)
    s=s.replace('**','').replace('__','').replace('`','')
    s=re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)', r'\1', s)
    return s.strip()

def annotate_inline_pmids(s):
    """Add PMID lookup labels directly beside numeric in-text citations."""
    def repl(m):
        raw=m.group(1)
        nums=re.findall(r'\d+', raw)
        if not nums:
            return m.group(0)
        bits=[]
        for n in nums:
            tag=pmids.get(int(n), UNVERIFIED)
            bits.append(f'{n} (PMID: {tag})')
        # Preserve the citation grouping while making every number searchable.
        return '[' + ', '.join(bits) + ']'
    return re.sub(r'\[([0-9][0-9,\-–\s]*)\]', repl, s)

def set_cell_text(cell, txt):
    cell.text=''
    p=cell.paragraphs[0]
    p.add_run(clean_inline(txt))

def set_repeat_table_header(row):
    trPr = row._tr.get_or_add_trPr(); tblHeader = OxmlElement('w:tblHeader'); tblHeader.set(qn('w:val'),'true'); trPr.append(tblHeader)

doc=Document()
sec=doc.sections[0]
sec.top_margin=Inches(0.8); sec.bottom_margin=Inches(0.8); sec.left_margin=Inches(0.9); sec.right_margin=Inches(0.9)
styles=doc.styles
styles['Normal'].font.name='Times New Roman'; styles['Normal']._element.rPr.rFonts.set(qn('w:eastAsia'),'Times New Roman'); styles['Normal'].font.size=Pt(10.5)
for sty in ['Heading 1','Heading 2','Heading 3']:
    styles[sty].font.name='Times New Roman'; styles[sty]._element.rPr.rFonts.set(qn('w:eastAsia'),'Times New Roman'); styles[sty].font.bold=True
styles['Heading 1'].font.size=Pt(14); styles['Heading 2'].font.size=Pt(12); styles['Heading 3'].font.size=Pt(11)

# Main body before References; omit explicit figure links, legend paragraphs, and
# Markdown tables.  Tables are supporting analysis and must not be inserted into
# this text-only manuscript Word file.
in_ref=False; skip_legend=False; table_buf=[]
def flush_table():
    # Deliberately discard Markdown table blocks.  They are not part of the
    # manuscript body and should be supplied separately as tables/results files.
    global table_buf
    table_buf=[]

for l in lines[:ref_start]:
    s=l.strip()
    if not s: flush_table(); continue
    if s.startswith('[Figure ') or s.startswith('![Figure'):
        flush_table(); continue
    if re.match(r'^(\*\*)?Fig\.\s*[0-9S]', s) or s.startswith('Figure S') or s=='Supplementary figure legends':
        flush_table(); continue
    if s.startswith('|'):
        table_buf.append(s); continue
    flush_table()
    if s.startswith('### '): doc.add_heading(clean_inline(s[4:]), level=3)
    elif s.startswith('## '): doc.add_heading(clean_inline(s[3:]), level=2)
    elif s.startswith('# '): doc.add_heading(clean_inline(s[2:]), level=1)
    else:
        body_text=annotate_inline_pmids(clean_inline(s))
        p=doc.add_paragraph(body_text); p.paragraph_format.space_after=Pt(5)
flush_table()

# References with PMID search labels.
doc.add_page_break(); doc.add_heading('References', level=1)
for n,body in refs:
    tag=pmids.get(n, UNVERIFIED)
    p=doc.add_paragraph(); p.paragraph_format.left_indent=Inches(0.15); p.paragraph_format.first_line_indent=Inches(-0.15); p.paragraph_format.space_after=Pt(4)
    # Put the PMID immediately after the reference number so it can be searched
    # manually in EndNote; it is part of each reference entry, not a separate
    # block placed before the reference list.
    p.add_run(f'{n}. [PMID: {tag}] ').bold=True
    p.add_run(clean_inline(body))

# Figure legends after references, without embedding figures.
doc.add_page_break(); doc.add_heading('Figure legends', level=1)
for lg in legends:
    p=doc.add_paragraph(clean_inline(lg)); p.paragraph_format.space_after=Pt(6)

# Footer with page number only.
for section in doc.sections:
    footer=section.footer.paragraphs[0]; footer.alignment=WD_ALIGN_PARAGRAPH.CENTER
    run=footer.add_run('Manuscript — '); fld=OxmlElement('w:fldSimple'); fld.set(qn('w:instr'),'PAGE'); footer._p.append(fld)

doc.core_properties.title='Immune cell composition shapes the interpretation of blood transcriptional recovery during tuberculosis treatment'
doc.core_properties.subject='Final text-only manuscript with PMID search labels for manual EndNote insertion'
doc.save(out)
print(json.dumps({'output':str(out),'references':len(refs),'legends':len(legends),'pmid_labels':len(pmids)},ensure_ascii=False))

