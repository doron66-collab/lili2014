#!/usr/bin/env python3
# Build a CGU-compliant PDF from dissertation_revised.html WITHOUT touching the source.
# CGU: 1" margins, 12pt Times, double-spaced body, roman prelim (title/committee/abstract
# counted-not-visible, first visible roman on TOC), arabic body from 1, page # bottom-center.
import re, subprocess, os, html as htmllib
from bs4 import BeautifulSoup
from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
import io

SRC   = '/home/user/lili2014/dissertation_revised.html'
OUT   = '/home/user/lili2014/dissertation_CGU.pdf'
SC    = '/tmp/claude-0/-home-user-lili2014/58416a4c-567e-5f83-9938-2ff40fbb2e74/scratchpad'
CHROME= '/opt/pw-browsers/chromium-1194/chrome-linux/chrome'
YEAR  = '2026'
AUTHOR= 'Doron Cohen'
TITLE = ('SOLANGE: A Provenance-Complete, Enterprise-Architected QC&middot;AI&middot;HPC '
         'Orchestration Platform for Non-Druggable NSCLC Tumor-Suppressor Mutations')

soup = BeautifulSoup(open(SRC, encoding='utf-8').read(), 'html.parser')

CSS = """
@page { size: letter; margin: 1in; }
* { box-sizing: border-box; }
body { font-family: 'Times New Roman', Times, serif; font-size: 12pt; line-height: 2;
       color: #000; margin: 0; }
h1 { font-size: 14pt; line-height: 2; page-break-before: always; margin: 0 0 12pt 0; }
h2 { font-size: 13pt; line-height: 2; margin: 12pt 0 4pt 0; }
h3 { font-size: 12pt; line-height: 2; margin: 10pt 0 4pt 0; }
p  { margin: 0; text-align: justify; text-indent: 0.4in; }
.center p, .titlepage p { text-indent: 0; }
.center { text-align: center; }
table { border-collapse: collapse; width: 100%; line-height: 1.15; font-size: 10pt;
        page-break-inside: avoid; margin: 8pt 0; }
th, td { border: 1px solid #000; padding: 3pt 5pt; text-align: left; vertical-align: top; }
th { font-weight: bold; }
ul { line-height: 2; margin: 0 0 0 0; }
.ref { line-height: 1.15; margin: 0 0 6pt 0; font-size: 12pt; padding-left: 0.3in; text-indent: -0.3in; }
.prelim { page-break-after: always; }
.titlepage { text-align: center; }
.spacer { height: 1.1in; }
.spacer-s { height: 0.5in; }
svg { max-width: 100%; height: auto; }
figure { page-break-inside: avoid; margin: 10pt 0; text-align: center; }
"""

def clean(s): return re.sub(r'\s+', ' ', s or '').strip()
def esc(s): return htmllib.escape(s)

def render(html_str, out_pdf):
    path = os.path.join(SC, 'frag.html')
    open(path, 'w', encoding='utf-8').write(
        f'<!doctype html><html><head><meta charset="utf-8"><style>{CSS}</style></head>'
        f'<body>{html_str}</body></html>')
    subprocess.run([CHROME, '--headless', '--no-sandbox', '--disable-gpu',
                    '--no-pdf-header-footer', f'--print-to-pdf={out_pdf}',
                    'file://' + path], capture_output=True, timeout=120)
    return out_pdf

# ── Build content fragments ──────────────────────────────────────────────────
def frag_title():
    return ('<div class="titlepage"><div class="spacer"></div>'
            f'<p><b>{TITLE}</b></p><div class="spacer-s"></div>'
            '<p>A dissertation submitted to the faculty of Claremont Graduate University<br>'
            'in partial fulfillment of the requirements for the degree of Doctor of Philosophy</p>'
            '<div class="spacer-s"></div><p>by</p>'
            f'<p><b>{AUTHOR}</b></p><div class="spacer-s"></div>'
            f'<p>Claremont Graduate University</p><p>{YEAR}</p></div>')

def frag_committee():
    return ('<div class="center"><p><b>Approval of the Dissertation Committee</b></p><br>'
            f'<p>This dissertation has been duly read, reviewed, and critiqued by the Committee '
            f'listed below, which hereby approves the manuscript of {AUTHOR} as fulfilling the '
            f'scope and quality requirements for meriting the degree of Doctor of Philosophy.</p>'
            '<br><br><p>Prof. Itamar Shabtai, Chair<br>Claremont Graduate University</p><br>'
            '<p>[Committee Member]</p><p>[Committee Member]</p></div>')

def frag_abstract():
    sec = soup.find('section', id='abstract')
    ps = ''.join(f'<p>{esc(clean(p.get_text(" ", strip=True)))}</p>' for p in sec.find_all('p')
                 if clean(p.get_text()))
    return ('<div class="center"><p><b>Abstract of the Dissertation</b></p>'
            f'<p><b>{TITLE}</b></p><p>by</p><p>{AUTHOR}</p>'
            f'<p>Claremont Graduate University: {YEAR}</p></div><br>' + ps)

SECTION_ORDER = ['problem','questions','literature','framework','methodology',
                 'contributions','preliminary','limitations','significance','timeline',
                 'references','appendix-a']

def frag_toc():
    rows = ''
    for sid in SECTION_ORDER:
        sec = soup.find('section', id=sid)
        if not sec: continue
        num = sec.find(class_='section-number'); title = sec.find(class_='section-title')
        n = clean(num.get_text(' ', strip=True)) if num else ''
        t = clean(title.get_text(' ', strip=True)) if title else ''
        rows += f'<p>{esc(n)} &mdash; {esc(t)}</p>'
    return ('<h1 style="page-break-before:avoid" class="center">Table of Contents</h1>'
            '<br>' + rows +
            '<br><p style="font-size:10pt"><i>(Page numbers to be finalized on conversion; '
            'this preview demonstrates CGU formatting &mdash; 1&Prime; margins, 12pt Times, '
            'double spacing.)</i></p>')

def render_table_html(tbl):
    out = '<table>'
    for r in tbl.find_all('tr'):
        cells = r.find_all(['th', 'td'])
        out += '<tr>'
        for c in cells:
            tag = 'th' if c.name == 'th' else 'td'
            out += f'<{tag}>{esc(clean(c.get_text(" ", strip=True)))}</{tag}>'
        out += '</tr>'
    return out + '</table>'

def frag_body():
    parts = []
    for sid in SECTION_ORDER:
        sec = soup.find('section', id=sid)
        if not sec: continue
        num = sec.find(class_='section-number'); title = sec.find(class_='section-title')
        n = clean(num.get_text(' ', strip=True)) if num else ''
        t = clean(title.get_text(' ', strip=True)) if title else ''
        parts.append(f'<h1>{esc(n)} &mdash; {esc(t)}</h1>')
        skip = set()
        for el in sec.find_all(['h3','h4','p','table','ul','ol','svg','div'], recursive=True):
            if any(a in skip for a in el.parents): continue
            cls = el.get('class', [])
            if el.name == 'div' and 'ref' in cls:
                parts.append(f'<p class="ref">{esc(clean(el.get_text(" ", strip=True)))}</p>')
                skip.add(el); continue
            if el.name in ('h3','h4'):
                lvl = 'h2' if el.name == 'h3' else 'h3'
                tx = clean(el.get_text(' ', strip=True))
                if tx: parts.append(f'<{lvl}>{esc(tx)}</{lvl}>')
            elif el.name == 'p':
                if el.find_parent('table') or el.find_parent(class_='ref'): continue
                tx = clean(el.get_text(' ', strip=True))
                if tx: parts.append(f'<p>{esc(tx)}</p>')
            elif el.name == 'table':
                parts.append(render_table_html(el)); skip.add(el)
            elif el.name in ('ul','ol'):
                if el.find_parent('table'): continue
                lis = ''.join(f'<li>{esc(clean(li.get_text(" ", strip=True)))}</li>'
                              for li in el.find_all('li', recursive=False) if clean(li.get_text()))
                if lis: parts.append(f'<ul>{lis}</ul>')
                skip.add(el)
            elif el.name == 'svg':
                parts.append(f'<figure>{str(el)}<br><span style="font-size:10pt">'
                             f'[Figure &mdash; rendered from source]</span></figure>')
                skip.add(el)
    return ''.join(parts)

# ── Render each block to its own PDF ─────────────────────────────────────────
blocks = {
    'title':     (frag_title(),     'roman', False),  # counted, not visible
    'committee': (frag_committee(), 'roman', False),
    'abstract':  (frag_abstract(),  'roman', False),
    'toc':       (frag_toc(),       'roman', True),    # first visible roman
    'body':      (frag_body(),      'arabic', True),
}
order = ['title','committee','abstract','toc','body']
pdfs = {}
for name in order:
    frag, _, _ = blocks[name]
    pdfs[name] = render(frag, os.path.join(SC, f'cgu_{name}.pdf'))

# ── Page-number overlay ──────────────────────────────────────────────────────
def roman(n):
    vals = [(1000,'m'),(900,'cm'),(500,'d'),(400,'cd'),(100,'c'),(90,'xc'),(50,'l'),
            (40,'xl'),(10,'x'),(9,'ix'),(5,'v'),(4,'iv'),(1,'i')]
    out=''
    for v,s in vals:
        while n>=v: out+=s; n-=v
    return out

def stamp(pdf_path, numbers):
    """numbers: list same length as pages; None = no number."""
    reader = PdfReader(pdf_path)
    writer = PdfWriter()
    for i, page in enumerate(reader.pages):
        num = numbers[i] if i < len(numbers) else None
        if num is not None:
            w = float(page.mediabox.width); h = float(page.mediabox.height)
            buf = io.BytesIO()
            c = canvas.Canvas(buf, pagesize=(w, h))
            c.setFont('Times-Roman', 11)
            c.drawCentredString(w/2, 0.5*72, str(num))  # 0.5in from bottom
            c.save(); buf.seek(0)
            overlay = PdfReader(buf).pages[0]
            page.merge_page(overlay)
        writer.add_page(page)
    tmp = pdf_path.replace('.pdf', '_num.pdf')
    with open(tmp, 'wb') as f: writer.write(f)
    return tmp

# assign numbers
counts = {name: len(PdfReader(pdfs[name]).pages) for name in order}
roman_counter = 1
final_parts = []
for name in order:
    n = counts[name]
    _, scheme, visible = blocks[name]
    if scheme == 'roman':
        nums = []
        for _ in range(n):
            nums.append(roman(roman_counter) if visible else None)
            roman_counter += 1
    else:  # arabic body, restart at 1
        nums = [str(i+1) for i in range(n)]
    final_parts.append(stamp(pdfs[name], nums))

# ── Merge ────────────────────────────────────────────────────────────────────
writer = PdfWriter()
for part in final_parts:
    for p in PdfReader(part).pages:
        writer.add_page(p)
with open(OUT, 'wb') as f: writer.write(f)
print('saved', OUT, '| total pages:', sum(counts.values()))
print('breakdown:', counts)
