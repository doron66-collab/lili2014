const fs = require('fs');
const { Document, Packer, Paragraph, TextRun, AlignmentType, Header, PageNumber, LineRuleType } = require('docx');

const F = 'Times New Roman';
const sp = { line: 480, lineRule: LineRuleType.AUTO, before: 0, after: 0 };

// Parse *italic* and **bold** markers into runs
function runs(text, base = {}) {
  const out = [];
  const re = /(\*\*[^*]+\*\*|\*[^*]+\*)/g;
  let last = 0, m;
  while ((m = re.exec(text))) {
    if (m.index > last) out.push(new TextRun({ text: text.slice(last, m.index), font: F, size: 24, ...base }));
    const t = m[0];
    if (t.startsWith('**')) out.push(new TextRun({ text: t.slice(2, -2), font: F, size: 24, bold: true, ...base }));
    else out.push(new TextRun({ text: t.slice(1, -1), font: F, size: 24, italics: true, ...base }));
    last = m.index + t.length;
  }
  if (last < text.length) out.push(new TextRun({ text: text.slice(last), font: F, size: 24, ...base }));
  return out;
}
const center = (t, o = {}) => new Paragraph({ alignment: AlignmentType.CENTER, spacing: sp, keepNext: !!o.keepNext, pageBreakBefore: !!o.pb, children: runs(t, o.bold ? { bold: true } : {}) });
const blank = () => new Paragraph({ spacing: sp, children: [new TextRun({ text: '', font: F, size: 24 })] });
const body = (t) => new Paragraph({ spacing: sp, indent: { firstLine: 720 }, children: runs(t) });
const h1 = (t) => new Paragraph({ alignment: AlignmentType.CENTER, spacing: sp, keepNext: true, children: runs(t, { bold: true }) });
const ref = (t) => new Paragraph({ spacing: sp, indent: { left: 720, hanging: 720 }, children: runs(t) });

const content0 = JSON.parse(fs.readFileSync('content.json', 'utf8'));
const title = content0.title;
const content = JSON.parse(fs.readFileSync('content.json', 'utf8'));

const titlePage = [
  blank(), blank(), blank(),
  center(title, { bold: true }),
  blank(),
  center('Doron Cohen'),
  center('Center for Information Systems and Technology, Claremont Graduate University'),
  center('TNDY 315: Principles of Project Management in a Complex World'),
  center('Dr. Matthew Muga'),
  center(content.date),
];

const bodyParas = [center(title, { bold: true, pb: true })];
for (const b of content.body) bodyParas.push(b.h ? h1(b.h) : body(b.p));

const refs = [center('References', { bold: true, pb: true })];
for (const r of content.refs) refs.push(ref(r));

const doc = new Document({
  styles: { default: { document: { run: { font: F, size: 24 } } } },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    headers: { default: new Header({ children: [new Paragraph({ alignment: AlignmentType.RIGHT, children: [new TextRun({ children: [PageNumber.CURRENT], font: F, size: 24 })] })] }) },
    children: [...titlePage, ...bodyParas, ...refs],
  }],
});
Packer.toBuffer(doc).then(b => fs.writeFileSync(content.out, b));
