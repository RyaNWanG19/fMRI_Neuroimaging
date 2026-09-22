import fs from 'node:fs/promises';
import path from 'node:path';
import {pathToFileURL} from 'node:url';
import {Presentation,PresentationFile} from '@oai/artifact-tool';
const dir=path.dirname(new URL(import.meta.url).pathname.replace(/^\/(\w:)/,'$1'));
const workspaceDir=path.dirname(dir), out=path.join(workspaceDir,'output');
const skill='C:/Users/Ryan_/.codex/plugins/cache/openai-primary-runtime/presentations/26.905.11957/skills/presentations';
const {finalizePresentation,applyPresentationChartFont}=await import(pathToFileURL(path.join(skill,'container_tools/artifact_tool_utils.mjs')));
const d=JSON.parse(await fs.readFile(path.join(dir,'results.json'),'utf8'));
// Chart workbooks use Excel's decimal precision. Retain source JSON at full
// precision and round plotted means to 9 decimal places, times to 2.
d.timeSec=d.timeSec.map(v=>Number(v.toFixed(2)));
d.examples.forEach(e=>e.mean=e.mean.map(row=>row.map(v=>Number(v.toFixed(9)))));
await fs.mkdir(out,{recursive:true});
const p=Presentation.create({slideSize:{width:1280,height:720}});
const font='Arial', ink='#192B36', muted='#53636D', colors=['#1675A9','#D47726','#6A4AA2'];
const chartFont={typeface:font,fontSize:22,fill:ink};
const source=d.models.map(m=>`${m.name}: ${m.source}`).join('\n');
const commonNotes=`Source: project results, Body_LoadDiff_vs_Face_LoadDiff.\n${source}\n${d.notes}\n99,999 sign flips and 99,999 centered bootstrap replicates. Holm alpha 0.05 separately across 489 ROIs within each model. Bootstrap is approximate. Model overlaps are descriptive, not study-wide multiplicity control.\n`;
function text(s,str,x,y,w,h,size=26,bold=false,color=ink){
 const t=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 t.text=str; t.text.style={typeface:font,fontSize:size,bold,color,autoFit:'none'}; return t;
}
function slide(title,subtitle,notes){const s=p.slides.add();s.background.fill='#FFFFFF';text(s,title,50,28,1180,66,40,true);if(subtitle)text(s,subtitle,52,100,1176,57,24,false,muted);s.speakerNotes.textFrame.setText(commonNotes+notes);return s;}
function foot(s,str){text(s,str,52,655,1176,46,19,false,muted);}
function chart(s,type,config){const c=s.charts.add(type,{chartFill:'#FFFFFF',plotAreaFill:'#FFFFFF',...config});applyPresentationChartFont(c,{fontFamily:font});return c;}
function table(s,values,x,y,width,height,widths){
 const t=s.tables.add({rows:values.length,columns:values[0].length,left:x,top:y,width,height,columnWidths:widths,values});
 t.borders.assign({fill:'#DDE4E8',width:1,style:'solid'});
 t.cells.block({row:0,column:0,rowCount:values.length,columnCount:values[0].length}).assign({textStyle:{typeface:font,fontSize:24,color:ink},fill:'#FFFFFF',margins:{left:12,right:12,top:10,bottom:10}});
 t.cells.block({row:0,column:0,rowCount:1,columnCount:values[0].length}).assign({fill:'#EAF1F5',textStyle:{typeface:font,fontSize:23,bold:true,color:ink}});
 return t;
}
// 1. Discovery counts and population interpretation.
{
const s=slide('Reliable residual contrasts across many ROIs','WM interaction: (2-back body − 0-back body) − (2-back faces − 0-back faces)',
'Talk track: Each bar counts regions whose population mean residual contrast differs from zero somewhere in the analyzed curve. cHRF detects 332 of 489, cHRFderiv 304 and sHRF 114. These are ROI-level results, not 41 independent lag tests. Discovery count differences do not constitute a test of HRF-model differences.');
chart(s,'bar',{position:{left:52,top:185,width:760,height:405},categories:d.models.map(m=>m.name),series:[{name:'Significant ROIs',values:d.models.map(m=>m.significant),points:colors.map((fill,idx)=>({idx,fill}))}],barOptions:{direction:'column',grouping:'clustered',gapWidth:100},hasLegend:false,xAxis:{textStyle:chartFont,majorGridlines:null},yAxis:{min:0,max:489,majorUnit:100,title:'ROIs passing Holm correction',textStyle:chartFont,majorGridlines:{fill:'#E1E6E9',width:1}},dataLabels:{showValue:true,position:'outEnd',textStyle:{...chartFont,bold:true,fontSize:27}}});
text(s,'67.9%   cHRF\n62.2%   cHRFderiv\n23.3%   sHRF',875,210,330,180,30,true);
text(s,'The population-average curve departs from zero in these regions.',875,415,315,145,27);
foot(s,'489 ROIs · 402–410 retained subjects per ROI · 41 lags · Holm p ≤ 0.05 within each model');
}
// 2. Exact disjoint overlap groups, avoiding misleading pairwise counts.
{
const s=slide('104 ROIs pass in all three HRF models','Each bar is a mutually exclusive pattern of significant models',
`Talk track: ${d.sharedAll} regions pass in all three models, while ${d.any} pass in at least one and ${d.none} pass in none. The largest additional group, 191 regions, passes in cHRF and cHRFderiv only. Overlap shows stability to model choice, not independent replication or a formal model comparison. Exact patterns ordered as all three, cHRF plus derivative only, cHRF plus sHRF only, derivative plus sHRF only, cHRF only, derivative only, sHRF only. Counts: ${d.patternCounts.join(', ')}.`);
chart(s,'bar',{position:{left:52,top:174,width:1174,height:425},categories:['All three','cHRF + deriv','cHRF + sHRF','deriv + sHRF','cHRF only','deriv only','sHRF only'],series:[{name:'ROIs',values:d.patternCounts,points:d.patternCounts.map((v,idx)=>({idx,fill:idx===0?'#1675A9':'#93A8B7'}))}],barOptions:{direction:'column',gapWidth:60},hasLegend:false,xAxis:{textStyle:{...chartFont,fontSize:19}},yAxis:{min:0,max:220,majorUnit:50,title:'Number of ROIs',textStyle:chartFont,majorGridlines:{fill:'#E1E6E9',width:1}},dataLabels:{showValue:true,position:'outEnd',textStyle:{...chartFont,bold:true}}});
foot(s,'343 ROIs pass in at least one model. 146 pass in none. “deriv” denotes cHRFderiv. Overlap is descriptive.');
}
// 3. Editable scatter curves with numerical time axis and zero reference.
{
const s=slide('Mean trajectories show the direction and shape','Illustrative ROIs selected after testing. Each curve is a population-average residual contrast.',
`Talk track: Bilateral V1 passes the omnibus test in every model. Right a24 passes only in cHRF. Mean trajectories retain the effect direction that the squared omnibus statistic cannot show. Different significant/not-significant outcomes do not establish a statistically significant difference between models. ${d.selectionNote} Y-axis units are those of the input residual contrast; the source metadata does not establish percent BOLD. Displayed times use the project's stored convention. No pointwise significance or uncertainty bands are shown. Atlas names use CANlab2018 ROI indices.`);
const names=['Left V1 (ROI 1)','Right V1 (ROI 2)','Right a24 (ROI 122)'];
for(let j=0;j<3;j++){
 const ex=d.examples[j],x=42+j*407;
 chart(s,'scatter',{position:{left:x,top:173,width:393,height:380},title:names[j],titleTextStyle:{...chartFont,fontSize:26,bold:true},series:[{name:'Zero',xValues:d.timeSec,values:d.timeSec.map(()=>0),line:{fill:'#8B9297',width:1},marker:{symbol:'none'}},...ex.mean.map((v,k)=>({name:d.models[k].name,xValues:d.timeSec,values:v,line:{fill:colors[k],width:2.5},marker:{symbol:'none'}}))],hasLegend:false,scatterOptions:{style:'line'},xAxis:{min:0,max:28.8,majorUnit:7.2,title:'Time (s)',textStyle:{...chartFont,fontSize:19},numberFormatCode:'0.0',majorGridlines:null},yAxis:{title:j===0?'Mean residual contrast':undefined,textStyle:{...chartFont,fontSize:18},numberFormatCode:'0.00',majorGridlines:{fill:'#E1E6E9',width:1}}});
 text(s,j<2?'Passes in all three models':'Passes in cHRF only',x+35,571,370,38,22,true);
}
colors.forEach((c,k)=>text(s,d.models[k].name,450+k*170,611,170,30,23,true,c));
foot(s,'Separate y-axis scales. Zero reference shown. Omnibus significance does not identify significant individual lags.');
}
// 4. Native table with exact agreement, no asymmetry proof claims.
{
const s=slide('Bootstrap sensitivity results closely agree','The centered bootstrap allows asymmetric subject curves and provides an approximate comparison.',
'Talk track: The centered bootstrap and sign flips agree on 486 of 489 ROI decisions in each canonical model and 488 of 489 in sHRF. Three cHRF and three derivative discoveries do not pass bootstrap; sHRF adds one bootstrap discovery. 103 ROIs pass both methods in all three models. This does not verify symmetry or independence. The primary analysis is sign flipping, rather than selecting the smaller p-value from the two methods. Both use 99,999 replicates and separately apply Holm across all 489 ROIs.');
table(s,[['HRF model','Sign flips','Bootstrap','Both methods','Agreement'],...d.models.map(m=>[m.name,''+m.significant,''+m.bootstrap,''+m.both,`${((489-m.disagreement)/489*100).toFixed(1)}%`])],64,194,1152,264,[250,215,215,235,237]);
text(s,'103 ROIs pass both methods in all three models',64,506,1150,57,31,true);
text(s,'Close agreement supports stability to the calibration method.',64,575,1150,42,26);
foot(s,'Counts refer to Holm p ≤ 0.05. Independence is owner-confirmed. Null symmetry remains an assumption.');
}
// 5. Native evidence table and clear limits.
{
const s=slide('ROI evidence and interpretation','Examples connect the statistical evidence to the trajectories shown above.',
'Talk track: Each p-value is adjusted within one 489-ROI model family. Many raw p-values reach the Monte Carlo floor of 0.00001, producing tied Holm p-values of 0.00489. The first two example ROIs therefore cannot be ranked precisely by their p-values. Right a24 is an example of model-dependent detection, not evidence of a significant between-model difference. The central conclusion is systematic population-average WM residual contrast. This test does not measure classifier importance or causal influence and does not establish significant lags. 12 cHRF, 9 derivative and 28 sHRF ROIs were flagged by the raw-p two-SE Monte Carlo sensitivity screen; that heuristic is not a confidence interval. Suggested method-slide title: Regions with reliable population-average residual contrasts. In the methods, D_i(t) should be the vector across ROIs; D_i,r(t) is its scalar ROI component. Atlas source: local load_atlas(canlab2018).');
table(s,[['ROI','cHRF p_Holm','cHRFderiv p_Holm','sHRF p_Holm'],...d.examples.map(e=>[`${e.name} (${e.id})`,...e.pHolm.map(v=>v.toFixed(5))])],64,188,1152,250,[342,250,310,250]);
text(s,'Supported conclusion',64,477,560,40,28,true);
text(s,'Systematic population-average residual contrasts persist across many ROIs.',64,529,545,95,26);
text(s,'Interpretation boundary',680,477,530,40,28,true);
text(s,'Timing, classifier importance and causal influence require separate analyses.',680,529,525,95,26);
foot(s,'Smallest raw p = 0.00001. Tied adjusted p-values limit fine ranking. Model differences require a direct test.');
}
const candidate=path.join(dir,'candidate.pptx');
await (await PresentationFile.exportPptx(p)).save(candidate);
for(let i=0;i<p.slides.items.length;i++){
 const s=p.slides.items[i], blob=await p.export({slide:s,format:'png',scale:1});
 await fs.writeFile(path.join(dir,`slide-${i+1}.png`),new Uint8Array(await blob.arrayBuffer()));
}
const final=path.join(out,'WM_residual_omnibus_results_final.pptx');
await finalizePresentation({workspaceDir,candidatePath:candidate,finalPath:final,
 pythonExecutable:'C:/Users/Ryan_/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe',
 integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit','--require-native-table-slide','4','--require-native-table-slide','5'],
 requiredNativeTableOwnerSlides:[4,5],requiredNativeChartOwnerSlides:[1,2,3],
 explicitTotalSlideCount:5,materializeLiteralChartWorkbooks:true,
 fontPolicy:{basis:'design',families:[font]},verifyArtifactToolImport:true,
 receiptPath:path.join(dir,'validation-final.json')});
console.log(final);
