const fs=require('fs');
const N=require('./deck_notes.js');
const {Document,Packer,Paragraph,TextRun,HeadingLevel,AlignmentType,PageBreak}=require('docx');
const F='Calibri';
const titles=['Title','The problem: cancer mutations with no drug','Why quantum computing, and why not yet','What SOLANGE does','SOLANGE today: a working platform','What earlier projects teach us','Our main model: Cynefin','Adding a communication model for people','Team roles and responsibilities','Project timeline','Plan walkthrough: where we are today','The hardest phase: 3B on IBM Heron r3','Phase 4: independent evaluation','Management tracks next to the research','Key project artifacts','What this course changed in our plan','Thank you'];
const kids=[new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:120},children:[new TextRun({text:'SOLANGE Final Presentation: Speaker Script by Slide',bold:true,size:32,font:F})]}),
 new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:300},children:[new TextRun({text:'Doron: slides 1-5 and 16-17  |  Manuel: slides 6-10  |  Amit: slides 11-15  |  Slide 18 (references) is not presented',size:20,font:F,italics:true})]})];
Object.keys(N).forEach((k,i)=>{
  const txt=N[k]; const lines=txt.split('\n'); const spk=lines[0];
  kids.push(new Paragraph({keepNext:true,spacing:{before:280,after:80},children:[new TextRun({text:`Slide ${i+1}: ${titles[i]}`,bold:true,size:26,font:F,color:'1B2A41'})]}));
  kids.push(new Paragraph({keepNext:true,spacing:{after:100},children:[new TextRun({text:spk,bold:true,size:21,font:F,color:'B8441F'})]}));
  lines.slice(1).join('\n').split(/\n\s*\n/).map(p=>p.trim()).filter(Boolean).forEach(p=>kids.push(new Paragraph({spacing:{after:120,line:300},children:[new TextRun({text:p,size:22,font:F})]})));
});
const doc=new Document({sections:[{properties:{page:{size:{width:12240,height:15840},margin:{top:1080,bottom:1080,left:1200,right:1200}}},children:kids}]});
Packer.toBuffer(doc).then(b=>{fs.writeFileSync('SOLANGE_Speaker_Script.docx',b);console.log('ok')});
