const fs=require('fs');
const {Document,Packer,Paragraph,TextRun,AlignmentType}=require('docx');
const data=JSON.parse(fs.readFileSync(process.argv[2],'utf8'));
const F='Calibri';
const kids=[new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:120},children:[new TextRun({text:'SOLANGE Final Presentation: Speaker Script by Slide',bold:true,size:32,font:F})]}),
 new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:80},children:[new TextRun({text:'Doron: slides 1-5 and 16-17  |  Manuel: slides 6-10  |  Amit: slides 11-15  |  Slide 18 (references) is not presented',size:20,font:F,italics:true})]}),
 new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:300},children:[new TextRun({text:'Generated directly from the speaker notes inside the PowerPoint file, so every slide number matches the deck.',size:18,font:F,italics:true,color:'6B7785'})]})];
data.forEach(({n,title,notes})=>{
  const lines=notes.split('\n'); const first=lines[0].trim();
  kids.push(new Paragraph({keepNext:true,spacing:{before:280,after:80},children:[new TextRun({text:`Slide ${n}: ${title}`,bold:true,size:26,font:F,color:'1B2A41'})]}));
  const body = first.startsWith('SPEAKER') ? lines.slice(1).join('\n') : notes;
  kids.push(new Paragraph({keepNext:true,spacing:{after:100},children:[new TextRun({text:first.startsWith('SPEAKER')?first:'NOT PRESENTED',bold:true,size:21,font:F,color:'B8441F'})]}));
  if (first.startsWith('SPEAKER')) body.split(/\n\s*\n/).map(p=>p.trim()).filter(Boolean).forEach(p=>kids.push(new Paragraph({spacing:{after:120,line:300},children:[new TextRun({text:p,size:22,font:F})]})));
});
const doc=new Document({sections:[{properties:{page:{size:{width:12240,height:15840},margin:{top:1080,bottom:1080,left:1200,right:1200}}},children:kids}]});
Packer.toBuffer(doc).then(b=>{fs.writeFileSync('SOLANGE_Speaker_Script.docx',b);console.log('ok')});
