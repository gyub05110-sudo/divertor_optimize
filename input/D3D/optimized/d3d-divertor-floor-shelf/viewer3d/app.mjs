import * as THREE from 'three';
import {OrbitControls} from './vendor/OrbitControls.js';
import {revolve,heatAtTriangle,heatSpan,densify,densifyWithHeat} from './geometry.mjs';
import {deviceOf,lowerBranchOf,chooseInitialCase,displaySettings,cameraSettings,formatChange,transportLabel,sourceDescription} from './device.mjs';
import {setupHistory} from './history-ui.mjs';
import {getQuality,attachQuality,showQuality,statusColors,statusLabel,patchAt} from './quality.mjs';

const $=id=>document.getElementById(id);
const params=new URLSearchParams(location.hash.slice(1));
const catalog=await fetch('./catalog.json').then(r=>{if(!r.ok)throw Error('평형 목록 로드 실패');return r.json();});
const initialEntry=chooseInitialCase(catalog.cases,params);
let data=await fetch('./'+initialEntry.bundle_url).then(r=>{if(!r.ok)throw Error('평형 데이터 로드 실패');return r.json();});
const initialBundle=data;
let qualityReport=await getQuality(initialEntry.id).catch(()=>null),qualityEpoch=0;
const state={model:'candidate',width:data.widths.includes(3)?3:data.widths[0],angle:240,wall:true,plasma:true,flux:true,sol:true,density:3};
let selectedEntry=initialEntry,stageMeta=data.stages.at(-1),historyController,unavailable=false;
if(['baseline','candidate'].includes(params.get('model')))state.model=params.get('model');
if(data.widths.includes(+params.get('width')))state.width=+params.get('width');
if(+params.get('angle')>=60&&+params.get('angle')<=360)state.angle=+params.get('angle');
$('width').value=state.width;$('angle').value=state.angle;
document.querySelector(`input[name=model][value=${state.model}]`).checked=true;
const palette=['#000004','#420a68','#932667','#dd513a','#fca50a','#fcffa4'];
const ramp=palette.map(c=>new THREE.Color(c));
function rgb(q,max){
 const t=Math.max(0,Math.min(1,q/max))*(ramp.length-1),i=Math.min(ramp.length-2,Math.floor(t));
 return ramp[i].clone().lerp(ramp[i+1],t-i);
}
const model=()=>data.models[state.model];
const maximum=()=>Math.max(data.models.baseline.summary[state.width].all_wall_peak_q_over_P,data.models.candidate.summary[state.width].all_wall_peak_q_over_P);
const auditedTarget=(m,s)=>(m.auditRanges||[]).some(([a,b])=>s>=a&&s<=b);
function sampledWall(m){return deviceOf(data)==='cpd2'||data.source_model==='unified_radial_exponential_v1'?densifyWithHeat(m.wall,m.heat[state.width]):densify(m.wall);}
let dense={baseline:sampledWall(data.models.baseline),candidate:sampledWall(data.models.candidate)};
const targets=m=>m.targets?m.targets.map(t=>t.points):[m.floor,m.shelf];
let renderer,scene,camera,controls,group,pickMeshes=[],pending=false,rendererFailure=null;
const host=$('scene');

function triangles(points,angle,knownOnly=false,sectors=88){
 const n=Math.max(12,Math.round(angle/360*sectors)),positions=[],colors=[],heatValues=[];
 const samples=model().heat[state.width],max=maximum();
 for(let i=1;i<points.length;i++){
  const a=points[i-1],b=points[i];
  if(knownOnly&&auditedTarget(model(),(a[2]+b[2])/2))continue;
  const span=knownOnly?heatSpan(samples,a[2],b[2]):[0,0];
  if(!span)continue;
  const [qa,qb]=span;
  for(let j=0;j<n;j++){
   const phi0=j/n*angle*Math.PI/180,phi1=(j+1)/n*angle*Math.PI/180;
   const corners=[revolve(a[0],a[1],phi0),revolve(b[0],b[1],phi0),revolve(b[0],b[1],phi1),revolve(a[0],a[1],phi1)];
   const qs=[qa,qb,qb,qa];
   for(const v of [0,1,2,0,2,3]){positions.push(...corners[v]);heatValues.push(qs[v]);if(knownOnly)colors.push(...rgb(qs[v],max).toArray());}
  }
 }
 const geo=new THREE.BufferGeometry();geo.setAttribute('position',new THREE.Float32BufferAttribute(positions,3));
 if(knownOnly){geo.setAttribute('color',new THREE.Float32BufferAttribute(colors,3));geo.setAttribute('q',new THREE.Float32BufferAttribute(heatValues,1));}
 geo.computeVertexNormals();return geo;
}
function outline(points,phi,color,opacity=1){
 const geo=new THREE.BufferGeometry().setFromPoints(points.map(p=>new THREE.Vector3(...revolve(p[0],p[1],phi))));
 return new THREE.Line(geo,new THREE.LineBasicMaterial({color,transparent:opacity<1,opacity}));
}
function auditedPatches(patches,angle){
 const n=Math.max(12,Math.round(angle/360*88)),positions=[],colors=[],heatValues=[];
 for(const p of patches){
  const q=p.q[state.width],color=q===null?new THREE.Color(statusColors[p.status]||'#ed65b3'):rgb(q,maximum());
  for(let j=0;j<n;j++){
   const f0=j/n*angle*Math.PI/180,f1=(j+1)/n*angle*Math.PI/180;
   const corners=[revolve(...p.a,f0),revolve(...p.b,f0),revolve(...p.b,f1),revolve(...p.a,f1)];
   for(const v of [0,1,2,0,2,3]){positions.push(...corners[v]);colors.push(...color.toArray());heatValues.push(q===null?NaN:q);}
  }
 }
 const geo=new THREE.BufferGeometry();geo.setAttribute('position',new THREE.Float32BufferAttribute(positions,3));geo.setAttribute('color',new THREE.Float32BufferAttribute(colors,3));geo.setAttribute('q',new THREE.Float32BufferAttribute(heatValues,1));
 return new THREE.Mesh(geo,new THREE.MeshBasicMaterial({vertexColors:true,side:THREE.DoubleSide,polygonOffset:true,polygonOffsetFactor:-3,polygonOffsetUnits:-3}));
}
function clearGroup(){if(!group)return;group.traverse(o=>{o.geometry?.dispose();if(o.material)for(const mat of Array.isArray(o.material)?o.material:[o.material])mat.dispose();});scene.remove(group);}
function rebuild(){
 clearGroup();group=new THREE.Group();pickMeshes=[];
 const pts=dense[state.model],angle=state.angle;
 if(state.wall){
  const mesh=new THREE.Mesh(triangles(pts,angle),new THREE.MeshBasicMaterial({color:0x92a0ae,transparent:true,opacity:.14,side:THREE.DoubleSide,depthWrite:false}));group.add(mesh);
  for(const phi of [0,angle*Math.PI/180])group.add(outline(model().wall,phi,0x8896a3,.85));
  for(const p of model().wall.filter((_,i)=>i%8===0)){
   const ring=[];for(let j=0;j<=64;j++)ring.push(new THREE.Vector3(...revolve(p[0],p[1],angle*Math.PI/180*j/64)));
   group.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(ring),new THREE.LineBasicMaterial({color:0xabb6c2,transparent:true,opacity:.25})));
  }
 }
 // Show the design targets even where first-hit samples do not cover them.
 for(const points of targets(model())){
  const target=densify(points);
  group.add(new THREE.Mesh(triangles(target,angle),new THREE.MeshBasicMaterial({color:0x748391,transparent:true,opacity:.55,side:THREE.DoubleSide,depthWrite:false})));
  group.add(outline(points,0,0x596978));group.add(outline(points,angle*Math.PI/180,0x596978));
 }
 const heat=new THREE.Mesh(triangles(pts,angle,true),new THREE.MeshBasicMaterial({vertexColors:true,side:THREE.DoubleSide,polygonOffset:true,polygonOffsetFactor:-2,polygonOffsetUnits:-2}));
 heat.userData.heat=true;group.add(heat);pickMeshes.push(heat);
 if(model().audit){const patches=auditedPatches(model().audit.patches,angle);group.add(patches);pickMeshes.push(patches);}
 if(state.plasma){
  const lcfs=[...data.lcfs,data.lcfs[0]];
  group.add(new THREE.Mesh(triangles(densify(lcfs,.02),angle),new THREE.MeshBasicMaterial({color:0x3586b0,transparent:true,opacity:.055,side:THREE.DoubleSide,depthWrite:false})));
  for(const phi of [0,angle*Math.PI/180])group.add(outline(lcfs,phi,0x3b87b2,.8));
 }
 if(state.flux){
  const surfaces=state.density===3?data.flux_surfaces.filter((_,i)=>[0,2,4].includes(i)):data.flux_surfaces;
  for(const surface of surfaces){
   const pts=thinCurve(surface.points,150);
   group.add(new THREE.Mesh(triangles(pts,angle,false,32),new THREE.MeshBasicMaterial({color:0x8d63bf,transparent:true,opacity:.027,side:THREE.DoubleSide,depthWrite:false})));
   for(const phi of [0,angle*Math.PI/180])group.add(outline(pts,phi,0x8d63bf,.65));
  }
 }
 if(state.sol){
  const paths=model().sol.filter((_,i)=>state.density===5||i%2===0);
  for(const path of paths){
   const pts=thinCurve(path.points,150);
   group.add(new THREE.Mesh(triangles(pts,angle,false,28),new THREE.MeshBasicMaterial({color:0x298764,transparent:true,opacity:.025,side:THREE.DoubleSide,depthWrite:false})));
   for(const phi of [0,angle*Math.PI/180])group.add(outline(pts,phi,0x288f67,.72));
  }
 }
 // A tiny orientation tripod, independent of the data heat scale.
 const axes=new THREE.AxesHelper(.6);axes.position.set(0,-1.55,0);group.add(axes);
 scene.add(group);invalidate();
}
function invalidate(){if(!renderer||pending)return;pending=true;requestAnimationFrame(()=>{pending=false;controls.update();renderer.render(scene,camera);});}
function resize(){if(renderer){const w=host.clientWidth,h=host.clientHeight;renderer.setSize(w,h,false);camera.aspect=w/h;camera.updateProjectionMatrix();invalidate();}draw2D();}
function preset(name){
 if(!camera)return;
 const c=cameraSettings(data);controls.maxDistance=c.maxDistance;
 if(name==='focus'){controls.target.set(...c.focus);camera.position.set(...c.focusPosition);}
 else if(name==='section'){controls.target.set(...c.section);camera.position.set(...c.sectionPosition);}
 else if(name==='top'){controls.target.set(0,c.z,0);camera.position.set(...c.top);}
 else{controls.target.set(0,c.z,0);camera.position.set(...c.overview);}
 controls.update();invalidate();
}
function zoom(factor){if(!camera)return;camera.position.sub(controls.target).multiplyScalar(factor).add(controls.target);controls.update();invalidate();}
function update(){
 if(unavailable)return;
 const b=data.models.baseline.summary[state.width],o=data.models.candidate.summary[state.width];
 const settings=displaySettings(data);
 if((deviceOf(data)==='d3d'||data.source_model==='unified_radial_exponential_v1')&&data.sigma_over_lambda_q!=null)$('widthNote').textContent=transportLabel(data)+` · 선택한 폭의 σ_OMP=${(state.width*data.sigma_over_lambda_q).toLocaleString('en-US')} mm · 확산 조건은 실행 목록에서 선택하세요.`;
 $('basePeak').textContent=b.all_wall_peak_q_over_P.toFixed(settings.peakDigits);
 $('optPeak').textContent=o.all_wall_peak_q_over_P.toFixed(settings.peakDigits);
 const change=100*(1-o.all_wall_peak_q_over_P/b.all_wall_peak_q_over_P);
 document.querySelector('.metrics .accent span').textContent=change>=0?'기준 대비 최대값 감소율':'기준 대비 최대값 증가율';
 $('reduction').textContent=formatChange(change);
 $('modelLabel').textContent=stageMeta.label;
 showQuality(qualityReport,model(),stageMeta,data);
 $('scaleMax').textContent=maximum().toFixed(3);$('angleValue').textContent=state.angle+'°';
 history.replaceState(null,'','#'+new URLSearchParams({device:deviceOf(data),case:selectedEntry.id,stage:stageMeta.id,model:state.model,width:state.width,angle:state.angle}));
 dense={baseline:sampledWall(data.models.baseline),candidate:sampledWall(data.models.candidate)};
 if(renderer)rebuild();draw2D();
 window.viewerState={...state,deviceId:deviceOf(data),units:settings.unit,caseId:selectedEntry.id,stageId:stageMeta.id,stageCount:data.stages.length,sigmaOverLambda:data.sigma_over_lambda_q??0,evaluationBudget:data.evaluation_budget??null,baselinePeak:b.all_wall_peak_q_over_P,candidatePeak:o.all_wall_peak_q_over_P,qualityApplied:!!model().audit,qualityVerified:qualityReport?.verified??false,qualitySourceCells:qualityReport?.levels.at(-1).source_cells,sampleCount:model().validation?.sample_count,heatVertexCount:pickMeshes[0]?.geometry?.getAttribute('q')?.count,ready:!!renderer};
}

function draw2D(){
 const canvas=$('sectionCanvas'),dpr=Math.min(devicePixelRatio||1,1.75),w=canvas.clientWidth,h=canvas.clientHeight;
 canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);
 const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);ctx.clearRect(0,0,w,h);
 if(unavailable)return;
 const stack=w<520,pw=stack?w:w/2,ph=stack?h/2:h;
 for(let index=0;index<2;index++){
  const key=index?'candidate':'baseline',m=data.models[key],x0=stack?0:index*pw,y0=stack?index*ph:0;
  const [rmin,rmax,zmin,zmax]=displaySettings(data).sectionBounds;
  const margin={l:40,r:15,t:35,b:28},scale=Math.min((pw-margin.l-margin.r)/(rmax-rmin),(ph-margin.t-margin.b)/(zmax-zmin));
  const X=r=>x0+margin.l+(r-rmin)*scale,Y=z=>y0+margin.t+(zmax-z)*scale;
  ctx.font='13px "Malgun Gothic",sans-serif';ctx.fillStyle='#34465b';ctx.fillText(index?'선택 단계 형상':'기준 형상',x0+margin.l,y0+19);
  ctx.save();ctx.beginPath();ctx.rect(x0+margin.l,y0+margin.t,pw-margin.l-margin.r,ph-margin.t-margin.b);ctx.clip();
  function line(points,color,lw,dash=[]){ctx.beginPath();points.forEach((p,i)=>i?ctx.lineTo(X(p[0]),Y(p[1])):ctx.moveTo(X(p[0]),Y(p[1])));ctx.strokeStyle=color;ctx.lineWidth=lw;ctx.lineCap='round';ctx.setLineDash(dash);ctx.stroke();ctx.setLineDash([]);}
  line(m.wall,'#b0bac4',1.5);
  if(state.plasma)line([...data.lcfs,data.lcfs[0]],'#3b87b2',1,[3,3]);
  if(state.flux)for(const c of data.flux_surfaces)line(c.points,'#9a73c5',.8);
  if(state.sol)for(const c of m.sol)line(c.points,'#369d74',.9);
  for(const points of targets(m))line(points,index?'#c77527':'#2166ac',2,index?[]:[5,3]);
  const pts=dense[key],ss=m.heat[state.width];
  for(let j=1;j<pts.length;j++){
   const a=pts[j-1],b=pts[j],span=heatSpan(ss,a[2],b[2]);
   if(auditedTarget(m,(a[2]+b[2])/2))continue;
   if(span)line([a,b],'#'+rgb((span[0]+span[1])/2,maximum()).getHexString(),5);
  }
  if(m.audit)for(const p of m.audit.patches){const q=p.q[state.width];line([p.a,p.b],q===null?statusColors[p.status]:'#'+rgb(q,maximum()).getHexString(),5);}
  const peak=m.summary[state.width];ctx.beginPath();ctx.arc(X(peak.peak_R_m),Y(peak.peak_Z_m),6,0,Math.PI*2);ctx.strokeStyle='#00a1b2';ctx.lineWidth=2;ctx.stroke();ctx.restore();
  ctx.fillStyle='#718096';ctx.font='11px sans-serif';ctx.fillText('R (m)',x0+pw-45,y0+ph-7);
  for(let i=0;i<4;i++){const r=rmin+(rmax-rmin)*(i+.2)/4;if(X(r)<x0+pw-10)ctx.fillText(r.toFixed(1),X(r)-8,y0+ph-15);}
 }
}

try{
 renderer=new THREE.WebGLRenderer({antialias:true,preserveDrawingBuffer:true});renderer.setPixelRatio(Math.min(devicePixelRatio||1,1.75));renderer.setClearColor(0xf8fafc);host.prepend(renderer.domElement);
 scene=new THREE.Scene();camera=new THREE.PerspectiveCamera(42,1,.03,200);
 controls=new OrbitControls(camera,renderer.domElement);controls.enableDamping=false;controls.minDistance=.35;controls.maxDistance=18;controls.addEventListener('change',invalidate);
 const ray=new THREE.Raycaster(),pointer=new THREE.Vector2();
 renderer.domElement.addEventListener('pointermove',ev=>{
  if(ev.buttons){$('tooltip').hidden=true;return;}
  const bounds=renderer.domElement.getBoundingClientRect();pointer.set((ev.clientX-bounds.left)/bounds.width*2-1,-(ev.clientY-bounds.top)/bounds.height*2+1);ray.setFromCamera(pointer,camera);
  const hit=ray.intersectObjects(pickMeshes)[0];
  if(!hit){$('tooltip').hidden=true;return;}
  const geo=hit.object.geometry,positions=geo.getAttribute('position'),values=geo.getAttribute('q'),ids=[hit.face.a,hit.face.b,hit.face.c];
  const local=hit.object.worldToLocal(hit.point.clone());
  const sample=heatAtTriangle(ids.map(i=>[positions.getX(i),positions.getY(i),positions.getZ(i)]),ids.map(i=>values.getX(i)),local.toArray());
  if(!sample){$('tooltip').hidden=true;return;}
  const {R:r,Z:z,q}=sample,patch=patchAt(model().audit?.patches,r,z);$('tooltip').hidden=false;$('tooltip').style.left=Math.max(0,Math.min(bounds.width-240,ev.clientX-bounds.left+14))+'px';$('tooltip').style.top=Math.max(55,ev.clientY-bounds.top-62)+'px';
  $('tooltip').textContent=`R ${r.toFixed(3)} m · Z ${z.toFixed(3)} m | ${patch?statusLabel(patch.status,patch.scope_reason)+' | ':''}q/P ${q===null?'미정':q.toPrecision(4)+' m⁻²'}`+(patch?.angle_deg!=null?` · 입사각 ${patch.angle_deg.toFixed(3)}°`:'');
 });
 renderer.domElement.addEventListener('pointerleave',()=>{$('tooltip').hidden=true;});
 $('status').hidden=true;preset('overview');resize();
}catch(error){renderer=null;rendererFailure='이 브라우저에서 3D를 초기화하지 못했습니다. 아래 2D 비교는 계속 확인할 수 있습니다.';$('status').textContent=rendererFailure;console.error(error);}

document.querySelectorAll('input[name=model]').forEach(el=>el.addEventListener('change',()=>historyController.selectModel(el.value)));
$('width').addEventListener('change',()=>{state.width=+$('width').value;update();});
$('angle').addEventListener('change',()=>{state.angle=+$('angle').value;update();});
$('angle').addEventListener('input',()=>{$('angleValue').textContent=$('angle').value+'°';});
for(const key of ['wall','plasma','flux','sol'])$(key).addEventListener('change',()=>{state[key]=$(key).checked;update();});
$('density').addEventListener('change',()=>{state.density=+$('density').value;update();});
for(const key of ['overview','focus','section','top'])$(key).addEventListener('click',()=>preset(key));
$('zoomIn').addEventListener('click',()=>zoom(.8));$('zoomOut').addEventListener('click',()=>zoom(1.25));
$('save').addEventListener('click',()=>{
 if(!renderer||unavailable)return;
 const filename=`${displaySettings(data).deviceLabel}-${selectedEntry.id}-${stageMeta.id}-${state.width}mm.png`;
 renderer.render(scene,camera);renderer.domElement.toBlob(blob=>{
  if(!blob)return;const url=URL.createObjectURL(blob),a=document.createElement('a');
  a.download=filename;a.href=url;a.hidden=true;document.body.append(a);a.click();a.remove();
  setTimeout(()=>URL.revokeObjectURL(url),10000);
 },'image/png');
});
host.addEventListener('keydown',ev=>{
 if(!camera)return;const delta=camera.position.clone().sub(controls.target);
 if(ev.key==='ArrowLeft'||ev.key==='ArrowRight')delta.applyAxisAngle(new THREE.Vector3(0,1,0),ev.key==='ArrowLeft'?.12:-.12);
 else if(ev.key==='ArrowUp'||ev.key==='ArrowDown')delta.y+=ev.key==='ArrowUp'?.25:-.25;
 else if(ev.key==='+'||ev.key==='='){zoom(.8);ev.preventDefault();return;}
 else if(ev.key==='-'){zoom(1.25);ev.preventDefault();return;}else return;
 camera.position.copy(controls.target).add(delta);controls.update();invalidate();ev.preventDefault();
});
function thinCurve(points,max){
 const step=Math.max(1,Math.ceil(points.length/max)),out=points.filter((_,i)=>i%step===0);
 if(out.at(-1)!==points.at(-1))out.push(points.at(-1));return out;
}
function applyStage(bundle,detail,entry,descriptor,viewParams){
 const token=++qualityEpoch;
 const changedDevice=deviceOf(data)!==deviceOf(bundle);
 if(qualityReport?.case_id!==bundle.id)qualityReport=null;
 const matching=qualityReport?.case_id===bundle.id&&qualityReport.shot_sha256===bundle.shot_sha256&&qualityReport.geometry_sha256===bundle.source_geometry_sha256?qualityReport:null;
 const base=attachQuality(bundle.models.baseline,matching?.models.baseline);
 const candidate=attachQuality(attachQuality(attachQuality(detail,detail.embedded_audit),matching?.models.candidate),matching?.models.baseline);
 data={...bundle,models:{baseline:base,candidate}};
 selectedEntry=entry;stageMeta=descriptor;state.model=descriptor.kind==='baseline'?'baseline':'candidate';unavailable=false;
 if(viewParams&&bundle.widths.includes(+viewParams.get('width')))state.width=+viewParams.get('width');
 if(!bundle.widths.includes(state.width))state.width=bundle.widths.includes(3)?3:bundle.widths[0];
 if(viewParams&&+viewParams.get('angle')>=60&&+viewParams.get('angle')<=360)state.angle=+viewParams.get('angle');
 $('angle').value=state.angle;
 $('width').replaceChildren();for(const width of bundle.widths){const option=document.createElement('option');option.value=width;option.textContent=Number(width).toFixed(Number.isInteger(width)?0:6);$('width').append(option);}
 $('width').value=state.width;$('width').disabled=bundle.widths.length===1;
 const settings=displaySettings(bundle);
 document.title=`${settings.deviceLabel} · 3D 열유속 비교`;
 document.querySelector('.eyebrow').textContent=`${settings.deviceLabel} / FIXED EQUILIBRIUM`;
 document.querySelector(`input[name=model][value=${state.model}]`).checked=true;
 document.querySelector('.tag').textContent=settings.tag;
 $('equilibriumFile').accept=settings.fileAccept;
 $('normalizationNote').textContent=sourceDescription(bundle);
 $('angleNote').textContent=settings.angleText;
 $('sol').parentElement.lastChild.textContent=lowerBranchOf(bundle)?'SOL 표시 (외측 하부 가지)':'SOL 표시 (양방향 폴로이달 가지)';
 $('solColorNote').textContent=lowerBranchOf(bundle)?'초록: 외측 하부 SOL · 표면 색: 열유속':'초록: 양방향 SOL · 표면 색: 열유속';
 $('legendUnits').textContent=settings.unit+' · 두 형상 공통 범위';
 $('solNote').textContent=lowerBranchOf(bundle)?'SOL은 외측 상류 중간면에서 하부 첫 벽 충돌까지의 폴로이달 경로를 회전한 면입니다.':'SOL은 저장된 양방향 경로를 실제 첫 벽 충돌에서 자른 폴로이달 경로의 회전면입니다.';
 $('widthNote').textContent=bundle.sigma_over_lambda_q!=null&&deviceOf(bundle)==='d3d'?transportLabel(bundle)+' · 확산 조건은 평형 / 실행 목록에서 선택하세요.':bundle.widths.length===1?'이 실행에서 계산한 λq입니다. 확산 가정은 평형 / 실행 목록에서 선택하세요.':'';
 $('fluxNote').textContent=deviceOf(bundle)==='cpd2'?'LCFS와 내부 자기면은 저장된 CPD 고정 평형에서 계산했습니다. ψ_N은 자기축 0 / 저장 경계 1입니다.':'LCFS는 저장된 EFIT 경계입니다. 내부 자기면은 ψ_N = 0.2·0.4·0.6·0.8·0.95의 닫힌 등고면입니다.';
 if(changedDevice)preset('overview');
 $('status').hidden=!!renderer;if(rendererFailure)$('status').textContent=rendererFailure;
 $('save').disabled=!renderer;update();
 if(deviceOf(bundle)==='d3d'&&qualityReport?.case_id!==bundle.id)getQuality(bundle.id).then(report=>{if(token!==qualityEpoch||!report)return;qualityReport=report;applyStage(bundle,detail,entry,descriptor);}).catch(error=>console.error(error));
}
function unavailableView(message){
 unavailable=true;clearGroup();group=null;pickMeshes=[];invalidate();draw2D();
 $('status').textContent=message;$('status').hidden=false;
 for(const id of ['basePeak','optPeak','reduction'])$(id).textContent='—';
 $('modelLabel').textContent='연결된 결과 없음';window.viewerState={ready:false,unavailable:true};
 qualityEpoch++;const panel=$('numericalQuality');if(panel)panel.textContent='연결된 수치검증 결과 없음';
}
new ResizeObserver(resize).observe(host);new ResizeObserver(draw2D).observe($('sectionCanvas'));
historyController=setupHistory({catalog,initialEntry,initialBundle,initialStage:params.get('stage')||(state.model==='baseline'?initialBundle.stages[0].id:initialBundle.stages.at(-1).id),applyStage,unavailableView});
await historyController.ready;
