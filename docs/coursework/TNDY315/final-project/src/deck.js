const pptxgen = require('pptxgenjs');
const React = require('react');
const ReactDOMServer = require('react-dom/server');
const sharp = require('sharp');
const fa = require('react-icons/fa');
const fs = require('fs');

const NAVY = '1B2A41', TERRA = 'B8441F', SLATE = '5A6B7E', LIGHT = 'F2F4F7', WHITE = 'FFFFFF',
  TEAL = '2E7D6B', GOLD = 'C28A1B', PURPLE = '6B4FA0', INK = '1F2933', MUTED = '6B7785';
const HF = 'Cambria', BF = 'Calibri';
const NOTES = require('./deck_notes.js');

async function icon(name, color, size = 256) {
  const svg = ReactDOMServer.renderToStaticMarkup(React.createElement(fa[name], { color: '#' + color, size }));
  const buf = await sharp(Buffer.from(svg)).png().toBuffer();
  return 'image/png;base64,' + buf.toString('base64');
}

(async () => {
  const pres = new pptxgen();
  pres.layout = 'LAYOUT_16x9';
  pres.title = 'SOLANGE: Planning a Quantum-Oncology Research Platform';
  pres.author = 'Doron Cohen, Manuel Blanco, Amit Mekel';

  const T = (s, text, o = {}) => s.addText(text, { isTextBox: true, fontFace: BF, color: INK, fontSize: 14, margin: 0, valign: 'top', ...o });
  const title = (s, text, sub) => {
    T(s, text, { x: 0.5, y: 0.32, w: 9, h: 0.6, fontFace: HF, fontSize: 28, bold: true, color: NAVY, valign: 'middle' });
    if (sub) T(s, sub, { x: 0.5, y: 0.9, w: 9, h: 0.35, fontSize: 13, italic: true, color: MUTED, valign: 'middle' });
  };
  let n = 0;
  const footer = (s, who, dark) => {
    n += 1;
    T(s, 'Presenter: ' + who, { x: 0.5, y: 5.22, w: 4, h: 0.25, fontSize: 9, color: dark ? 'AAB6C4' : MUTED, valign: 'middle' });
    T(s, String(n), { x: 9.0, y: 5.22, w: 0.5, h: 0.25, fontSize: 9, align: 'right', color: dark ? 'AAB6C4' : MUTED, valign: 'middle' });
  };
  const circleIcon = async (s, name, x, y, d, bg, fg = WHITE) => {
    s.addShape(pres.shapes.OVAL, { x, y, w: d, h: d, fill: { color: bg }, line: { color: bg } });
    s.addImage({ data: await icon(name, fg), x: x + d * 0.24, y: y + d * 0.24, w: d * 0.52, h: d * 0.52 });
  };
  const card = (s, x, y, w, h, fill = LIGHT) => s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, fill: { color: fill }, line: { color: fill }, rectRadius: 0.08 });
  const note = (s, key) => s.addNotes(NOTES[key]);

  // 1 Title
  let s = pres.addSlide(); s.background = { color: NAVY };
  await circleIcon(s, 'FaAtom', 0.6, 0.7, 0.8, TERRA);
  T(s, 'SOLANGE', { x: 0.6, y: 1.65, w: 8.8, h: 0.8, fontFace: HF, fontSize: 44, bold: true, color: WHITE });
  T(s, 'Planning a transdisciplinary quantum-oncology research platform under technological uncertainty', { x: 0.6, y: 2.45, w: 8.2, h: 0.9, fontSize: 20, color: 'D9E1EA' });
  T(s, 'Doron Cohen  |  Manuel Blanco  |  Amit Mekel', { x: 0.6, y: 3.65, w: 8.8, h: 0.35, fontSize: 15, bold: true, color: WHITE });
  T(s, 'TNDY 315: Principles of Project Management in a Complex World  |  Dr. Matthew Muga  |  Claremont Graduate University  |  September 2026', { x: 0.6, y: 4.05, w: 8.8, h: 0.5, fontSize: 11, color: 'AAB6C4' });
  footer(s, 'Doron Cohen', true); note(s, 's1');

  // 2 Problem
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'The problem: cancer mutations with no drug', 'Lung cancer is the most diagnosed and most lethal cancer worldwide');
  const stats = [['2.5M', 'new lung cancer cases in 2022', 'Bray et al., 2024'], ['1.8M', 'lung cancer deaths in 2022', 'Bray et al., 2024'], ['~60%', 'of NSCLC patients carry mutations with no approved targeted drug', 'Cohen, 2026a (estimate)']];
  stats.forEach(([big, lab, src], i) => {
    const x = 0.5 + i * 3.05;
    card(s, x, 1.5, 2.85, 2.1);
    T(s, big, { x: x + 0.2, y: 1.65, w: 2.45, h: 0.9, fontFace: HF, fontSize: 44, bold: true, color: i === 2 ? TERRA : NAVY, valign: 'middle' });
    T(s, lab, { x: x + 0.2, y: 2.6, w: 2.45, h: 0.65, fontSize: 13, color: INK });
    T(s, src, { x: x + 0.2, y: 3.25, w: 2.45, h: 0.25, fontSize: 9, italic: true, color: MUTED });
  });
  await circleIcon(s, 'FaMicroscope', 0.5, 3.95, 0.55, NAVY);
  T(s, 'Why: designing a drug needs an accurate model of the electrons at the mutated site. Classical computers lose accuracy when many electrons interact strongly.', { x: 1.25, y: 3.95, w: 8.25, h: 0.7, fontSize: 14, valign: 'middle' });
  footer(s, 'Doron Cohen'); note(s, 's2');

  // 3 Why quantum, why not yet
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'Why quantum computing, and why not yet', 'The technology we need matures after our milestones are due');
  // roadmap line
  s.addShape(pres.shapes.LINE, { x: 1.4, y: 2.2, w: 7.2, h: 0, line: { color: SLATE, width: 2 } });
  const road = [[1.4, '2023', 'IBM "utility" claim, later matched classically'], [3.8, '2026', 'Today: noisy hardware, pipeline tests on 3 IBM devices'], [6.2, '2027', 'Our defense (Dec 2027)'], [8.6, '2029', 'IBM Starling: first large-scale fault-tolerant system']];
  road.forEach(([x, yr, txt], i) => {
    s.addShape(pres.shapes.OVAL, { x: x - 0.14, y: 2.06, w: 0.28, h: 0.28, fill: { color: i === 3 ? TERRA : NAVY }, line: { color: WHITE, width: 2 } });
    T(s, yr, { x: x - 0.7, y: 1.55, w: 1.4, h: 0.4, fontFace: HF, fontSize: 18, bold: true, color: i === 3 ? TERRA : NAVY, align: 'center', valign: 'middle' });
    T(s, txt, { x: x - 0.85, y: 2.45, w: 1.7, h: 0.9, fontSize: 11, align: 'center', color: INK });
  });
  card(s, 0.5, 3.55, 9, 1.35);
  await circleIcon(s, 'FaExclamationTriangle', 0.75, 3.8, 0.8, TERRA);
  T(s, 'A wicked problem', { x: 1.8, y: 3.7, w: 7.5, h: 0.35, fontSize: 16, bold: true, color: NAVY });
  T(s, 'We cannot fully define what the platform must compute until we try to compute it (Rittel & Webber, 1972). Each discipline also defines success differently.', { x: 1.8, y: 4.05, w: 7.5, h: 0.75, fontSize: 13 });
  footer(s, 'Doron Cohen'); note(s, 's3');

  // 4 What SOLANGE does
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'What SOLANGE does', 'It decides, per mutation, whether a quantum computer is needed at all');
  const flow = [['FaDna', 'Mutation', 'from real patient sequencing data', NAVY], ['FaServer', 'Two classical methods', 'DMRG and SHCI run independently', SLATE], ['FaBalanceScale', 'Verdict', 'Class B: classically tractable\nClass A: quantum-necessary', TEAL], ['FaLock', 'LEON seals the record', 'P1 to P9 provenance, re-verified (FDA Part 11 logic)', TERRA]];
  for (let i = 0; i < flow.length; i++) {
    const [ic, h, b, col] = flow[i]; const x = 0.5 + i * 2.3;
    card(s, x, 1.5, 2.05, 2.35);
    await circleIcon(s, ic, x + 0.72, 1.65, 0.6, col);
    T(s, h, { x: x + 0.12, y: 2.35, w: 1.81, h: 0.5, fontSize: 14, bold: true, color: NAVY, align: 'center', valign: 'middle' });
    T(s, b, { x: x + 0.12, y: 2.85, w: 1.81, h: 0.9, fontSize: 11, align: 'center' });
    if (i < 3) s.addShape(pres.shapes.LINE, { x: x + 2.06, y: 2.67, w: 0.22, h: 0, line: { color: SLATE, width: 1.5, endArrowType: 'triangle' } });
  }
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.5, y: 4.15, w: 9, h: 0.8, fill: { color: 'FBEDE7' }, line: { color: 'FBEDE7' }, rectRadius: 0.08 });
  T(s, 'Finding so far: every real biological target we measured was classically tractable. The platform reports this openly.', { x: 0.75, y: 4.15, w: 8.5, h: 0.8, fontSize: 14, bold: true, color: TERRA, valign: 'middle' });
  footer(s, 'Doron Cohen'); note(s, 's4');

  // 4b Platform today
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'SOLANGE today: a working platform', 'Screens from the live system, September 2026');
  const pw = 5.2, ph = pw * 706 / 1351;
  s.addImage({ path: 'plat_plane.png', x: 0.5, y: 1.4, w: pw, h: ph });
  T(s, 'Classifier map: every real target so far lands in Class B (classically tractable)', { x: 0.5, y: 1.4 + ph + 0.08, w: pw, h: 0.45, fontSize: 11, italic: true, color: MUTED });
  const rw = 3.6;
  const h1 = rw * 490 / 1135, h2 = rw * 480 / 1178;
  s.addImage({ path: 'plat_qpu_s.png', x: 5.9, y: 1.4, w: rw, h: h1 });
  T(s, 'Test jobs on real IBM Heron hardware', { x: 5.9, y: 1.4 + h1 + 0.04, w: rw, h: 0.25, fontSize: 11, italic: true, color: MUTED });
  const y2 = 1.4 + h1 + 0.45;
  s.addImage({ path: 'plat_compliance_s.png', x: 5.9, y: y2, w: rw, h: h2 });
  T(s, 'LEON audit trail, P1 to P9 provenance', { x: 5.9, y: y2 + h2 + 0.04, w: rw, h: 0.25, fontSize: 11, italic: true, color: MUTED });
  footer(s, 'Doron Cohen'); note(s, 's4b');

  // 5 Prior work
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'What earlier projects teach us', 'Progress came when evidence set the pace, not marketing');
  const cols = [
    ['FaTimesCircle', 'What went wrong', TERRA, [['IBM Watson for Oncology', 'MD Anderson spent $62M, then canceled. Little evidence of patient benefit (Strickland, 2019).'], ['Quantum advantage claims', "IBM's 2023 result was later matched on classical computers, partly on a phone (Schirber, 2024)."]]],
    ['FaCheckCircle', 'What went well', TEAL, [['KRAS G12C', 'Pocket found in 2013, drug approved in 2021: eight years of patient work (Ostrem et al., 2013; FDA, 2021).'], ['Long partnerships', 'Cleveland Clinic and IBM: a 10-year program that also trains people (Cleveland Clinic, 2023).']]],
  ];
  for (let c = 0; c < 2; c++) {
    const [ic, h, col, items] = cols[c]; const x = 0.5 + c * 4.6;
    await circleIcon(s, ic, x, 1.45, 0.5, col);
    T(s, h, { x: x + 0.65, y: 1.45, w: 3.6, h: 0.5, fontSize: 17, bold: true, color: col, valign: 'middle' });
    items.forEach(([a, b], k) => {
      const y = 2.1 + k * 1.45;
      card(s, x, y, 4.4, 1.3);
      T(s, a, { x: x + 0.2, y: y + 0.12, w: 4.0, h: 0.35, fontSize: 14, bold: true, color: NAVY });
      T(s, b, { x: x + 0.2, y: y + 0.48, w: 4.0, h: 0.75, fontSize: 12 });
    });
  }
  footer(s, 'Manuel Blanco'); note(s, 's5');

  // 6 Cynefin
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'Our main model: Cynefin', 'PMBOK Guide, Section 4.2.5.1: match each part of the problem to its context');
  const q = [['Complicated', 'Known unknowns. Analyze, then apply good practice.', 'Classification stage: predictive plan, fixed milestones', NAVY, 0.5, 1.4],
    ['Complex', 'Unknown unknowns. Probe, sense, respond in short cycles.', 'Quantum stage: short cycles, quarterly re-plan', TERRA, 5.05, 1.4],
    ['Clear', 'Obvious cause and effect. Use best practice.', 'Routine operations', SLATE, 0.5, 3.2],
    ['Chaotic', 'Act first to stabilize, then sense.', 'Not expected in this project', SLATE, 5.05, 3.2]];
  q.forEach(([h, d, sol, col, x, y]) => {
    card(s, x, y, 4.45, 1.65, col === SLATE ? 'F7F8FA' : LIGHT);
    T(s, h, { x: x + 0.2, y: y + 0.12, w: 4.05, h: 0.4, fontSize: 17, bold: true, color: col });
    T(s, d, { x: x + 0.2, y: y + 0.52, w: 4.05, h: 0.5, fontSize: 12 });
    T(s, 'SOLANGE: ' + sol, { x: x + 0.2, y: y + 1.05, w: 4.05, h: 0.5, fontSize: 12, bold: true, color: col === SLATE ? MUTED : col });
  });
  footer(s, 'Manuel Blanco'); note(s, 's6');

  // 7 Communication model
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'Adding a communication model for people', 'Cynefin does not cover stakeholders (PMBOK Table 4-1)');
  const sh = [['FaFlask', 'Chemist', 'Converged energy'], ['FaUserMd', 'Oncologist', 'Clinical relevance'], ['FaCoins', 'Investor', 'Credible progress'], ['FaGlobeAmericas', 'Health ministry', 'Trust and access']];
  for (let i = 0; i < sh.length; i++) {
    const [ic, h, b] = sh[i]; const x = 0.5 + i * 2.3;
    await circleIcon(s, ic, x + 0.65, 1.45, 0.75, [NAVY, TEAL, GOLD, PURPLE][i]);
    T(s, h, { x, y: 2.3, w: 2.05, h: 0.35, fontSize: 14, bold: true, color: NAVY, align: 'center' });
    T(s, 'Success = ' + b, { x, y: 2.65, w: 2.05, h: 0.35, fontSize: 11, color: MUTED, align: 'center' });
  }
  card(s, 0.5, 3.25, 9, 1.7);
  await circleIcon(s, 'FaComments', 0.75, 3.55, 0.7, TERRA);
  T(s, 'Cross-cultural communication (PMBOK 4.2.2.1)', { x: 1.7, y: 3.38, w: 7.6, h: 0.35, fontSize: 15, bold: true, color: NAVY });
  T(s, [
    { text: 'Every result states its context: established or experimental.', options: { bullet: true, breakLine: true } },
    { text: 'Uncertain results are explained in a live briefing, not only in a report.', options: { bullet: true, breakLine: true } },
    { text: 'The same result must mean the same thing to every reader.', options: { bullet: true } },
  ], { x: 1.7, y: 3.75, w: 7.6, h: 1.1, fontSize: 13, paraSpaceAfter: 4 });
  footer(s, 'Manuel Blanco'); note(s, 's7');

  // 8 Roles
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'Team roles and responsibilities', 'Each member leads the area of his own professional expertise');
  const roles = [['FaAtom', 'Doron Cohen', 'Project lead and technical lead', ['Platform architecture and science', 'Research phases 1 to 7', 'Evidence discipline (LEON)'], NAVY],
    ['FaMoneyBillWave', 'Manuel Blanco', 'Funding and investor relations', ['Grants and binational funding', 'IBM Research Israel collaboration', 'Budget and funder reporting'], GOLD],
    ['FaHandshake', 'Amit Mekel', 'International partnerships', ['Health diplomacy', 'Partner hospitals and ministries', 'Expert-panel sourcing'], PURPLE]];
  for (let i = 0; i < 3; i++) {
    const [ic, nm, rl, items, col] = roles[i]; const x = 0.5 + i * 3.05;
    card(s, x, 1.45, 2.85, 3.5);
    await circleIcon(s, ic, x + 1.05, 1.65, 0.75, col);
    T(s, nm, { x: x + 0.15, y: 2.5, w: 2.55, h: 0.4, fontSize: 16, bold: true, color: NAVY, align: 'center' });
    T(s, rl, { x: x + 0.15, y: 2.9, w: 2.55, h: 0.35, fontSize: 12, italic: true, color: col, align: 'center' });
    T(s, items.map((t, k) => ({ text: t, options: { bullet: true, breakLine: k < items.length - 1 } })), { x: x + 0.25, y: 3.4, w: 2.4, h: 1.4, fontSize: 12, paraSpaceAfter: 4 });
  }
  footer(s, 'Manuel Blanco'); note(s, 's8');

  // 9 Timeline
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'Project timeline: March 2026 to February 2028');
  s.addImage({ path: 'gantt_crop.png', x: 1.166, y: 0.98, w: 7.669, h: 4.200 });
  footer(s, 'Manuel Blanco'); note(s, 's9');

  // 10 Walkthrough done
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'Plan walkthrough: where we are today', 'Phases 1 to 3A');
  const done = [['FaCheckCircle', 'Phase 1', 'Foundation and literature review', 'Complete', TEAL], ['FaCheckCircle', 'Phase 2', 'Working prototype deployed', 'Complete', TEAL], ['FaHourglassHalf', 'Phase 3A', 'Test the design on a simulator (May to Aug 2026)', 'Main run done Aug 4; larger run pending', TERRA], ['FaServer', 'Jul-Aug 2026', 'Small test jobs on 3 real IBM computers: Kingston, Marrakesh, Fez', 'Pipeline proven; not yet medical-grade chemistry', NAVY]];
  for (let i = 0; i < done.length; i++) {
    const [ic, h, d, st, col] = done[i]; const y = 1.45 + i * 0.9;
    card(s, 0.5, y, 9, 0.78);
    await circleIcon(s, ic, 0.65, y + 0.12, 0.54, col);
    T(s, h, { x: 1.35, y: y + 0.05, w: 1.35, h: 0.68, fontSize: 14, bold: true, color: NAVY, valign: 'middle' });
    T(s, d, { x: 2.7, y: y + 0.05, w: 3.9, h: 0.68, fontSize: 13, valign: 'middle' });
    T(s, st, { x: 6.7, y: y + 0.05, w: 2.65, h: 0.68, fontSize: 12, bold: true, color: col, valign: 'middle' });
  }
  footer(s, 'Amit Mekel'); note(s, 's10');

  // 11 Hardest phase
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'The hardest phase: 3B on IBM Heron r3', 'September to December 2026, contingent on access');
  const hp = [['Why it is hard', TERRA, 'FaExclamationTriangle', ['Access depends on IBM: 3 to 6 months', 'Immature hardware can fail on its own', 'Each discipline defines success differently']],
    ['How the plan responds', TEAL, 'FaRoute', ['Officially contingent: no access by December means future work, not failure', 'A target goes to hardware only if it fits what we already ran', 'Success criteria are written before each run']]];
  for (let c = 0; c < 2; c++) {
    const [h, col, ic, items] = hp[c]; const x = 0.5 + c * 4.6;
    card(s, x, 1.45, 4.4, 3.5);
    await circleIcon(s, ic, x + 0.25, 1.62, 0.55, col);
    T(s, h, { x: x + 0.95, y: 1.62, w: 3.3, h: 0.55, fontSize: 17, bold: true, color: col, valign: 'middle' });
    T(s, items.map((t, k) => ({ text: t, options: { bullet: true, breakLine: k < items.length - 1 } })), { x: x + 0.3, y: 2.4, w: 3.9, h: 2.4, fontSize: 14, paraSpaceAfter: 10 });
  }
  footer(s, 'Amit Mekel'); note(s, 's11');

  // 12 Phase 4
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'Phase 4: independent evaluation', 'January to April 2027, criteria frozen in advance on the Open Science Framework');
  const p4 = [['4A', 'Technical test harness', '4 wk'], ['4B', 'Comparative benchmarks', '3 wk'], ['4C', 'Recruit 4 to 6 experts', '3 wk'], ['4D', 'Two rounds of panel sessions', '4 wk'], ['4E', 'Real case walkthroughs', '2 wk'], ['4F', 'Triangulation and report', '3 wk']];
  p4.forEach(([id, d, wk], i) => {
    const x = 0.5 + i * 1.52;
    const col = id === '4C' ? PURPLE : NAVY;
    s.addShape(pres.shapes.OVAL, { x: x + 0.38, y: 1.55, w: 0.62, h: 0.62, fill: { color: col }, line: { color: col } });
    T(s, id, { x: x + 0.38, y: 1.55, w: 0.62, h: 0.62, fontSize: 15, bold: true, color: WHITE, align: 'center', valign: 'middle' });
    if (i < 5) s.addShape(pres.shapes.LINE, { x: x + 1.05, y: 1.86, w: 0.42, h: 0, line: { color: SLATE, width: 1.5, endArrowType: 'triangle' } });
    T(s, d, { x: x, y: 2.3, w: 1.38, h: 0.7, fontSize: 12, align: 'center', bold: true, color: NAVY });
    T(s, wk, { x: x, y: 3.0, w: 1.38, h: 0.3, fontSize: 11, align: 'center', color: MUTED });
  });
  card(s, 0.5, 3.6, 9, 1.35);
  await circleIcon(s, 'FaUsers', 0.75, 3.85, 0.8, PURPLE);
  T(s, 'The hardest step here is 4C: finding experts from several fields in three weeks. Our international partner network supports this recruitment.', { x: 1.8, y: 3.7, w: 7.45, h: 1.15, fontSize: 14, valign: 'middle' });
  footer(s, 'Amit Mekel'); note(s, 's12');

  // 13 Management tracks
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'Management tracks next to the research', 'People, money, and partnerships run in parallel with the science');
  const tr = [['FaSyncAlt', 'Quarterly re-plan', 'Every quarter, the plan is updated from what we actually achieved, never from vendor promises.', 'Whole team', TEAL],
    ['FaCoins', 'Funding track', 'Grants, binational funding, IBM Research Israel collaboration. October 2026 to June 2027.', 'Manuel', GOLD],
    ['FaGlobeAmericas', 'Health diplomacy', 'Partner hospitals and health ministries. Cancer does not stop at borders; shared science brings countries together.', 'Amit', PURPLE]];
  for (let i = 0; i < 3; i++) {
    const [ic, h, d, who, col] = tr[i]; const y = 1.45 + i * 1.2;
    card(s, 0.5, y, 9, 1.05);
    await circleIcon(s, ic, 0.7, y + 0.2, 0.65, col);
    T(s, h, { x: 1.6, y: y + 0.1, w: 2.3, h: 0.85, fontSize: 16, bold: true, color: col, valign: 'middle' });
    T(s, d, { x: 3.9, y: y + 0.1, w: 4.3, h: 0.85, fontSize: 12.5, valign: 'middle' });
    T(s, who, { x: 8.2, y: y + 0.1, w: 1.15, h: 0.85, fontSize: 12, bold: true, color: NAVY, align: 'right', valign: 'middle' });
  }
  footer(s, 'Amit Mekel'); note(s, 's13');

  // 14 Artifacts
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'Key project artifacts', 'Categories from PMBOK Guide Section 4.6');
  const ar = [['FaFileSignature', 'Project charter', 'Authority, milestones, 3B contingency'], ['FaMapSigns', 'Roadmap', 'Re-planned each quarter'], ['FaClipboardList', 'Assumption log', 'Every modeling choice behind a result'], ['FaExclamationTriangle', 'Risk register', 'Access, hardware, funding, recruitment'],
    ['FaUsers', 'Stakeholder register', 'What success means to each group'], ['FaExchangeAlt', 'Change log', 'Formal changes and retractions'], ['FaComments', 'Communications plan', 'Who hears what, and how'], ['FaMoneyBillWave', 'Budget', 'Compared with actuals each quarter']];
  for (let i = 0; i < ar.length; i++) {
    const [ic, h, d] = ar[i]; const x = 0.5 + (i % 4) * 2.3, y = 1.45 + Math.floor(i / 4) * 1.75;
    card(s, x, y, 2.1, 1.6);
    await circleIcon(s, ic, x + 0.15, y + 0.15, 0.5, [NAVY, NAVY, TEAL, TERRA, PURPLE, SLATE, PURPLE, GOLD][i]);
    T(s, h, { x: x + 0.15, y: y + 0.72, w: 1.85, h: 0.35, fontSize: 13, bold: true, color: NAVY });
    T(s, d, { x: x + 0.15, y: y + 1.05, w: 1.85, h: 0.5, fontSize: 11 });
  }
  footer(s, 'Amit Mekel'); note(s, 's14');

  // 15 Takeaways
  s = pres.addSlide(); s.background = { color: NAVY };
  T(s, 'What this course changed in our plan', { x: 0.5, y: 0.35, w: 9, h: 0.6, fontFace: HF, fontSize: 28, bold: true, color: WHITE, valign: 'middle' });
  const tk = [['FaRoute', 'Tailoring', 'One project, different approaches for different parts.'], ['FaBullseye', 'Models for coverage', 'Cynefin for uncertainty, communication for people.'], ['FaLightbulb', 'Evidence first', '"No quantum computer needed" is a result, not a failure.']];
  for (let i = 0; i < 3; i++) {
    const [ic, h, d] = tk[i]; const x = 0.5 + i * 3.05;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: 1.35, w: 2.85, h: 2.6, fill: { color: '26395A' }, line: { color: '26395A' }, rectRadius: 0.08 });
    await circleIcon(s, ic, x + 1.05, 1.55, 0.75, TERRA);
    T(s, h, { x: x + 0.2, y: 2.45, w: 2.45, h: 0.45, fontSize: 17, bold: true, color: WHITE, align: 'center' });
    T(s, d, { x: x + 0.2, y: 2.9, w: 2.45, h: 0.9, fontSize: 13, color: 'D9E1EA', align: 'center' });
  }
  T(s, 'Uncertain hardware becomes a managed condition: every phase delivers value on its own.', { x: 0.5, y: 4.3, w: 9, h: 0.6, fontSize: 15, italic: true, color: 'F3C9B8', align: 'center', valign: 'middle' });
  footer(s, 'Doron Cohen', true); note(s, 's15');

  // 16 Thank you
  s = pres.addSlide(); s.background = { color: NAVY };
  await circleIcon(s, 'FaAtom', 4.55, 0.9, 0.9, TERRA);
  T(s, 'Thank you', { x: 0.5, y: 2.0, w: 9, h: 0.8, fontFace: HF, fontSize: 40, bold: true, color: WHITE, align: 'center', valign: 'middle' });
  T(s, 'Questions and discussion', { x: 0.5, y: 2.8, w: 9, h: 0.45, fontSize: 18, color: 'D9E1EA', align: 'center' });
  T(s, 'Doron Cohen  |  Manuel Blanco  |  Amit Mekel', { x: 0.5, y: 3.6, w: 9, h: 0.4, fontSize: 14, bold: true, color: WHITE, align: 'center' });
  footer(s, 'Doron Cohen', true); note(s, 's16');

  // 17 References
  s = pres.addSlide(); s.background = { color: WHITE };
  title(s, 'References');
  const refs = require('./paper_content.js').refs.map(r => r.replace(/\*/g, ''));
  T(s, refs.map((r, k) => ({ text: r, options: { breakLine: k < refs.length - 1 } })), { x: 0.5, y: 1.0, w: 9, h: 4.15, fontSize: 7.5, paraSpaceAfter: 2, color: INK });
  footer(s, 'Not presented'); s.addNotes('Reference slide. Not presented aloud.');

  await pres.writeFile({ fileName: 'TNDY315_SOLANGE_Final_Presentation.pptx' });
  console.log('deck ok');
})();
