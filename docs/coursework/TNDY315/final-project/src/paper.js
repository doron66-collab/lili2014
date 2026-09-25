const fs = require('fs');
const C = require('./paper_content.js');
const { Document, Packer, Paragraph, TextRun, AlignmentType, Header, PageNumber, LineRuleType,
  PageOrientation, ImageRun, Table, TableRow, TableCell, WidthType, BorderStyle, PageBreak, VerticalAlign } = require('docx');

const F = 'Times New Roman';
const sp = { line: 480, lineRule: LineRuleType.AUTO, before: 0, after: 0 };
const sp1 = { line: 240, lineRule: LineRuleType.AUTO, before: 0, after: 0 };

function runs(text, base = {}, size = 24) {
  const out = []; const re = /(\*\*[^*]+\*\*|\*[^*]+\*)/g; let last = 0, m;
  while ((m = re.exec(text))) {
    if (m.index > last) out.push(new TextRun({ text: text.slice(last, m.index), font: F, size, ...base }));
    const t = m[0];
    if (t.startsWith('**')) out.push(new TextRun({ text: t.slice(2, -2), font: F, size, bold: true, ...base }));
    else out.push(new TextRun({ text: t.slice(1, -1), font: F, size, italics: true, ...base }));
    last = m.index + t.length;
  }
  if (last < text.length) out.push(new TextRun({ text: text.slice(last), font: F, size, ...base }));
  return out;
}
const P = (t, o = {}) => new Paragraph({ spacing: sp, indent: o.noIndent ? undefined : { firstLine: 720 }, alignment: o.align, children: runs(t, o.base || {}) });
const center = (t, o = {}) => new Paragraph({ alignment: AlignmentType.CENTER, spacing: sp, keepNext: true, pageBreakBefore: !!o.pb, children: runs(t, o.bold ? { bold: true } : {}) });
const blank = () => new Paragraph({ spacing: sp, children: [new TextRun({ text: '', font: F, size: 24 })] });
const h1 = (t, pb) => new Paragraph({ alignment: AlignmentType.CENTER, spacing: sp, keepNext: true, pageBreakBefore: !!pb, children: runs(t, { bold: true }) });
const h2 = (t) => new Paragraph({ spacing: sp, keepNext: true, children: runs(t, { bold: true }) });
const ref = (t) => new Paragraph({ spacing: sp, indent: { left: 720, hanging: 720 }, children: runs(t) });

const header = new Header({ children: [new Paragraph({ alignment: AlignmentType.RIGHT, children: [new TextRun({ children: [PageNumber.CURRENT], font: F, size: 24 })] })] });
const portrait = { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } };
const landscape = { page: { size: { width: 12240, height: 15840, orientation: PageOrientation.LANDSCAPE }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } };

// ---------- Section A: title page and body ----------
const A = [
  blank(), blank(), blank(),
  center(C.title, { bold: true }),
  blank(),
  center(C.authors),
  center('Center for Information Systems and Technology, Claremont Graduate University'),
  center('TNDY 315: Principles of Project Management in a Complex World'),
  center('Dr. Matthew Muga'),
  center(C.date),
  center(C.title, { bold: true, pb: true }),
  ...C.intro.map(t => P(t)),
  h1('Project Area: Quantum Computing for Non-Druggable Cancer Mutations'),
  ...C.s1.map(t => P(t)),
  h1('Prior Work in This Space'),
  h2('What Went Wrong'),
  ...C.s2a.map(t => P(t)),
  h2('What Went Well and Why'),
  ...C.s2b.map(t => P(t)),
  h1('How the Project Plan Was Constructed'),
  ...C.s3.map(t => P(t)),
  h2('Phases, Milestones, and Team Roles'),
  ...C.s3phases.map(t => P(t)),
  h2('The Most Difficult Phase'),
  ...C.s3hard.map(t => P(t)),
];

// ---------- Section B: landscape timeline ----------
const img = fs.readFileSync('gantt.png');
const B = [
  new Paragraph({ spacing: sp1, children: runs('Figure 1', { bold: true }) }),
  new Paragraph({ spacing: { ...sp1, after: 80 }, children: runs('*SOLANGE Project Timeline, March 2026 to February 2028*') }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: sp1, children: [new ImageRun({ data: img, type: 'png', transformation: { width: 864, height: Math.round(864 * 1639 / 2700) } })] }),
  new Paragraph({ spacing: { ...sp1, before: 60 }, children: runs('*Note.* Phases, durations, sub-phases, and milestones follow the dissertation timeline (Cohen, 2026a), updated to September 25, 2026. The management tracks and owners are the team plan for this project. Dashed bars are pending or contingent work. Phase 3B depends on IBM Quantum Network access.', {}, 20) }),
];

// ---------- Section C: artifacts table, model, references, appendix ----------
const W = [2050, 1450, 4760, 1100];
const line = { style: BorderStyle.SINGLE, size: 6, color: '000000' };
const none = { style: BorderStyle.NONE, size: 0, color: 'FFFFFF' };
function cell(text, w, o = {}) {
  return new TableCell({
    width: { size: w, type: WidthType.DXA },
    margins: { top: 50, bottom: 50, left: 80, right: 80 },
    verticalAlign: VerticalAlign.TOP,
    borders: { top: o.top ? line : none, bottom: o.bottom ? line : { style: BorderStyle.SINGLE, size: 2, color: 'BFBFBF' }, left: none, right: none },
    children: [new Paragraph({ spacing: sp1, children: runs(text, o.bold ? { bold: true } : {}, 20) })],
  });
}
const rowsT = [new TableRow({ tableHeader: true, children: ['Artifact', 'PMBOK category', 'How we would use it', 'Owner'].map((h, i) => cell(h, W[i], { bold: true, top: true, bottom: true })) })];
C.artifacts.forEach((r, k) => rowsT.push(new TableRow({ cantSplit: true, children: r.map((v, i) => cell(v, W[i], { bottom: k === C.artifacts.length - 1 })) })));
const table = new Table({ width: { size: 9360, type: WidthType.DXA }, columnWidths: W, rows: rowsT });

const shots = [
  ['plat_plane.png', 'SOLANGE platform, Dashboard: the classical tractability plane. Each point is a target classified by measured entanglement and active-space size (Cohen, 2026a).'],
  ['plat_classifier_c.png', 'SOLANGE platform, Orchestration: the two-method classifier (DMRG and SHCI) with sealed results for real targets (Cohen, 2026a).'],
  ['plat_qpu_c.png', 'SOLANGE platform, Orchestration: pipeline test jobs on real IBM Heron hardware, July to August 2026 (Cohen, 2026a).'],
  ['plat_compliance_c.png', 'SOLANGE platform, Compliance: the LEON audit trail and the P1 to P9 provenance schema (Cohen, 2026a).'],
  ['poster_hi-1.png', 'SOLANGE research poster, IST 697, May 2026.'],
  ['ms_p1-01.png', 'Manuscript prepared for arXiv submission, September 2026 working draft, first page (Cohen, 2026b).'],
  ['shot_title.png', 'SOLANGE dissertation proposal, title page (Cohen, 2026a).'],
  ['shot_abstract.png', 'SOLANGE dissertation abstract, stating that every real biological target measured to date is classically tractable (Cohen, 2026a).'],
];
const { imageSize } = (() => { try { return require('image-size'); } catch (e) { return {}; } })();
function dims(f, maxW, maxH) {
  const b = fs.readFileSync(f); const w = b.readUInt32BE(16), h = b.readUInt32BE(20);
  let s = Math.min(maxW / w, maxH / h); return { width: Math.round(w * s), height: Math.round(h * s) };
}
const appendix = [h1('Appendix: Supporting Artifacts', true), new Paragraph({ spacing: sp, children: runs('The screenshots below show the working SOLANGE platform and the dissertation behind the project plan.') })];
shots.forEach(([f, cap], i) => {
  appendix.push(new Paragraph({ spacing: { ...sp1, before: 200 }, keepNext: true, children: runs(`**Screenshot ${i + 1}.** ${cap}`, {}, 22) }));
  appendix.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: sp1, children: [new ImageRun({ data: fs.readFileSync(f), type: 'png', transformation: dims(f, 620, 400) })] }));
});

const Csec = [
  h1('Project Artifacts'),
  P(C.s4intro),
  new Paragraph({ spacing: sp1, keepNext: true, children: runs('Table 1', { bold: true }) }),
  new Paragraph({ spacing: { ...sp1, after: 80 }, keepNext: true, children: runs('*Artifacts the SOLANGE Project Would Maintain*') }),
  table,
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { ...sp, before: 240 }, keepNext: true, children: runs('Project Management Model', { bold: true }) }),
  ...C.s5.map(t => P(t)),
  h1('References', true),
  ...C.refs.map(ref),
  ...appendix,
];

const doc = new Document({
  styles: { default: { document: { run: { font: F, size: 24 } } } },
  sections: [
    { properties: portrait, headers: { default: header }, children: A },
    { properties: landscape, headers: { default: header }, children: B },
    { properties: portrait, headers: { default: header }, children: Csec },
  ],
});
Packer.toBuffer(doc).then(b => { fs.writeFileSync('TNDY315_SOLANGE_Final_Paper.docx', b); console.log('ok'); });
