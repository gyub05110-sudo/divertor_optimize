import * as THREE from 'three';
import {OrbitControls} from './vendor/OrbitControls.js';
import {revolve,segmentsForWall,selectedRecord,heatRGB,minRadius,maxTangentError,legLayers} from './view-math.mjs';
import {nearestHeatSegment,profileForView,segmentTriangleRZ,formatHeatQ,validateStageIndex,validateSavedIndex} from './heat-probe.mjs';
import {applyQualifiedFinals,validateQualifiedProfile} from './qualified-finals.mjs';
import {applyReassessment,distributionPolicy} from './distribution-checks.mjs';

const $=id=>document.getElementById(id);
const state={source:'IMP',radius:20,run:'campaign',index:0,final:false,zoom:false,profile:null};
let data,stages=null,savedProfiles=null,legs=null,legError='',series,record,timer=null,token=0,renderer,camera,controls,scene,group;
const legColors={inner:'#2270bd',outer:'#9d4387'};
const visibleLegs=()=> $('legs').checked?legLayers(legs,record,$('display').value):[];
const heatVisible=()=>!!state.profile&&$('heat').checked;
const heatMax=()=>data.baseline.refined[state.source];
const shownWall=()=>record&&$('display').value!=='baseline'?record.wall:data.baseline.wall;
let sectionFrame=null,pickWall=null,probeFrame=0,probePointer=null,profileLoading=false;
const probeRay=new THREE.Raycaster();
const profileCache=new Map();
const digestText=async text=>Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',new TextEncoder().encode(text))),v=>v.toString(16).padStart(2,'0')).join('');
const fixed=n=>Number(n).toFixed(4),percent=n=>(100*(1-n)).toFixed(2);
const distributionFailed=r=>r?.distribution_reassessment?.passed===false;
const finalLabel=r=>distributionFailed(r)?'이전 기준 통과 · v3 수치 기준 미통과':r?.distribution_policy===distributionPolicy?'v3 수치 기준 통과':r?.additional_validation?'추가 정밀검증 통과':r?.validation==='current'?'최종 검증 통과':'과거 정밀 결과';
const accent=()=>state.source==='IMP'?'#087e84':'#aa6b1b';
const error=e=>{$('error').hidden=false;$('error').textContent='기록을 표시하지 못했습니다. '+e.message;console.error(e);};

function stop(){clearInterval(timer);timer=null;$('play').textContent='▶';$('play').setAttribute('aria-label','개선 단계 자동 재생');}
function selectSeries(preferred=state.run){
  stop();
  const choices=data.series.filter(s=>s.source===state.source&&s.rmin===state.radius);
  series=choices.find(s=>s.run===preferred)??choices[0];
  state.run=series.run;state.index=series.stages.length;state.final=!!series.final;
  if(series.final?.additional_validation)state.index=series.final.archived_stage;
  $('runs').replaceChildren(...choices.map(s=>new Option(s.label,s.run)));
  $('runs').value=series.run;
  $('stage-select').replaceChildren(new Option('기준 형상','0'),...series.stages.map((r,i)=>new Option(`${String(i+1).padStart(2,'0')} · ${percent(r.ratio)}% 감소`,String(i+1))));
  $('timeline').max=series.stages.length;
  const stats=series.statistics;
  $('run-summary').textContent=`${stats.lines.toLocaleString()}개 평가 기록 · 표시 ${series.stages.length}단계${stats.excluded_geometry?` · 형상 검사 탈락 ${stats.excluded_geometry}단계`:''}`;
  $('final').disabled=!series.final;
  $('final').querySelector('span').textContent=series.final?`${finalLabel(series.final)} · ${percent(series.final.ratio)}% 감소`:'표시할 개선 정밀 결과 없음';
  document.body.dataset.source=state.source;
  for(const b of document.querySelectorAll('[data-source]'))if(b.tagName==='BUTTON')b.setAttribute('aria-pressed',String(b.dataset.source===state.source));
  for(const b of document.querySelectorAll('[data-radius]'))b.setAttribute('aria-pressed',String(Number(b.dataset.radius)===state.radius));
  update().catch(error);
}

async function update(){
  hideProbe();
  const version=++token;
  record=selectedRecord(series,state.index,state.final);state.profile=null;
  let rel=profileForView(data.baseline,state.source,record,$('display').value,stages);
  const saved=savedProfiles?.profiles?.[rel];
  const entry=saved??stages?.cases?.[record?.id];
  if(!entry?.sha256)throw new Error('해당 형상의 열분포 색인 누락');
  const expectedIdentity=saved?saved.identity:record.id,expectedDonors=saved?1280:320;
  if(saved)rel=saved.profile;
  profileLoading=!!rel;
  const base=record?.baseline??data.baseline.refined[state.source];
  $('case-label').textContent=`${state.source} / RMIN ${state.radius} MM`;
  $('stage-title').textContent=state.final?(distributionFailed(record)?'이전 최종 형상 · 국소 분포 차이 미해결':record.distribution_policy===distributionPolicy?'v3 수치 기준 통과 형상':record.additional_validation?'추가 정밀검증 결과':record.validation==='current'?'최종 수치검증 결과':'과거 정밀 결과'):record?`개선 단계 ${String(state.index).padStart(2,'0')}`:'D3D 기준 형상';
  $('base-peak').textContent=fixed(base);$('peak').textContent=fixed(record?.peak??base);
  $('reduction').textContent=record?`−${percent(record.ratio)}%`:'0.00%';
  $('resolution').textContent=`${state.final||!record?'1280':'320'} 경로 · 동일 해상도 기준 비교`;
  $('status').textContent=state.final?finalLabel(record):record?'탐색 기록 · 현재 형상 통과':'고정 기준';
  $('status').classList.toggle('verified',state.final&&record.validation==='current');
  $('timeline').value=state.index;$('stage-select').value=state.index;
  $('stage-counter').textContent=state.final?(record.additional_validation?`저장 ${record.archived_stage}단계 · 추가 정밀 평가`:'별도 정밀 평가'):`${state.index} / ${series.stages.length} 단계`;
  $('prev').disabled=state.index===0&&!state.final;
  $('next').disabled=state.index>=series.stages.length&&(!state.final||record?.additional_validation);
  $('play').disabled=series.stages.length===0;$('timeline').disabled=series.stages.length===0;
  $('empty').hidden=series.stages.length>0;
  $('empty').textContent='이 실행에는 현재 형상 검사를 통과한, 기준보다 개선된 기록 갱신 단계가 없습니다. 기준 형상을 표시합니다.';
  $('chart').hidden=series.stages.length===0;
  $('section-caption').textContent=state.zoom?`${state.source} 접합부 확대 · 왜곡 없는 동일 축척`:'동일 축척 · 기준은 점선, 선택 형상은 실선';
  $('section-zoom').textContent=state.zoom?'전체 단면':'접합부 확대';
  $('identity-short').textContent=record?record.id.slice(0,12):'original D3D';
  populateEvidence();drawChart();drawSection();buildScene();heatLegend();updateLegNote();
  if(rel){
    if(!profileCache.has(rel)){
      profileCache.set(rel,fetch(rel).then(async r=>{
        if(!r.ok)throw new Error('열유속 기록 읽기 실패');
        const text=await r.text();
        if(await digestText(text)!==entry.sha256)throw new Error('열분포 파일 해시 불일치');
        const p=JSON.parse(text);
        if(p.identity!==expectedIdentity||p.source!==entry.source||p.donors!==expectedDonors)throw new Error('열분포 후보 또는 해상도 불일치');
        if(entry.fingerprint)validateQualifiedProfile(p,entry);
        return p;
      }).catch(e=>{profileCache.delete(rel);throw e;}));
      while(profileCache.size>12)profileCache.delete(profileCache.keys().next().value);
    }
    let p;
    try {p=await profileCache.get(rel);}catch(e){if(version!==token)return;profileLoading=false;heatLegend();throw e;}
    if(version!==token)return;
    state.profile=p;profileLoading=false;heatLegend();drawSection();buildScene();
  }
}

function updateLegNote(){
  $('leg-legend').hidden=!$('legs').checked||!legs;
  const missing=record&&$('display').value!=='baseline'&&!legs?.cases[record.id];
  $('leg-note').textContent=legError||(!$('legs').checked?'Leg 표시 숨김':missing?'이 후보의 leg 기록이 없어 표시하지 않습니다.':
    'X-point → 각 벽의 첫 도착점 · 분리면의 R–Z 가지와 축대칭 회전면 · 나선형 이동 궤적이 아닙니다.');
}

function heatLegend(){
  const visible=heatVisible();
  document.querySelector('.heat-strip').classList.toggle('active',visible);
  const kind=record&&$('display').value!=='baseline'?'선택 형상':'기준 형상';
  const resolution=state.profile?.donors??1280;
  $('heat-note').textContent=state.profile?`${kind} · ${resolution} 경로${resolution===320?' · 해당 단계 재평가':''} · 표시 셀 평균${visible?' · ◇: peak 위치':' · 색상 숨김 / 마우스로 q 확인 가능'}`:profileLoading?'해당 형상의 열분포 불러오는 중…':'해당 형상의 열분포를 불러오지 못했습니다.';
  $('heat-scale').textContent=visible?'MW/m²':'';
  $('q-colorbar').hidden=!visible;
  const gradient=Array.from({length:21},(_,i)=>`rgb(${heatRGB(i/20,1).map(v=>Math.round(v*255)).join(',')}) ${i*5}%`).join(',');
  $('q-gradient').style.background=`linear-gradient(90deg,${gradient})`;
  $('q-ticks').replaceChildren(...[0,.25,.5,.75,1].map(t=>{const s=document.createElement('span');s.textContent=(heatMax()*t).toFixed(2);return s;}));
}

function populateEvidence(){
  const fields=record?[
    ['후보 ID',record.id],['평가 기록',state.final?'별도 정밀 결과':`순번 ${record.sequence} · 원본 ${record.line}행 · ${record.method}`],
    ['peak 해상도',record.resolution],['현재 형상 검사','접선·곡률·길이·전체 벽 교차 검사 통과'],
    ['실제 최소 반경',`${(1000*minRadius(record.constraints)).toFixed(2)} mm`],
    ['최대 접선 오차',`${maxTangentError(record.constraints).toExponential(2)}°`],
    ['열수송 검증',state.final?finalLabel(record):record.distribution_reassessment?(record.distribution_reassessment.passed?'이 후보의 분포 수렴 v3 검증 통과 · 화면 q는 320 경로':'이 후보의 분포 수렴 v3 미통과 · 탐색 기록 보존'):'이 단계의 정밀 열수송 검증은 확인되지 않음'],
    ['저장 열분포',record.profile?'해당 형상과 일치하는 저장 NPZ 사용':stages?.cases?.[record.id]?'동일 입력·320 경로로 이 후보 재평가 · 원래 peak 재현 확인':'열분포 없음'],
    ['원본 파일',record.origin]
  ]:[['형상','수정하지 않은 D3D 기준 벽'],['입력','IMP/OMP 각각 독립 1 MW · 고정 수송 입력'],['peak 해상도','cubic · 1280 donor'],['출처','campaign_geometry_20261029.json / 저장된 기준 정밀 평가']];
  if(record?.additional_validation)fields.push(
    ['추가 검증 시각',new Date(record.qualified_utc).toLocaleString('ko-KR',{timeZone:'Asia/Seoul'})+' (한국시간)'],
    ['원래 탐색 단계',`${record.archived_stage} / ${series.stages.length}단계 · 새 형상 탐색 없이 재검증`],
    ['선택 범위',distributionFailed(record)?'이전 검증 정책의 형상 보존 · peak 감소와 국소 분포 신뢰도는 별도 판단':'보관 후보의 1280 경로 peak 순위 중 해당 기준 통과 후보 · 미세한 순위 차이의 우월성은 미확정'],
    ['보관 기록 출처',record.archived_origin]);
  if(record?.distribution_reassessment)fields.push(
    ['분포 수렴 v3',record.distribution_reassessment.passed?'설정된 수치 기준 통과':'국소 q 차이 기준 미통과 · 개선 없음이나 물리적 부적합을 뜻하지 않음'],
    ['위치별 q 차이',`${(record.distribution_reassessment.max_linf_relative*100).toFixed(3)}% · 전역 peak 대비 · 잠정 수치 허용오차 1%`],
    ['허용오차의 근거','1%·0.5%는 연구용 수치 설정 · 물리적 허용조건에서 도출하거나 오차 상한으로 인증한 값 아님'],
    ['최대점 위치','진단값으로만 기록 · 봉우리 순위 교체만으로 탈락시키지 않음']);
  const dl=document.createElement('dl');
  for(const [key,val] of fields){const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=key;dd.textContent=val;dl.append(dt,dd);}
  $('evidence-body').replaceChildren(dl);
}

function drawChart(){
  const ratios=[1,...series.stages.map(s=>s.ratio)],maxDrop=Math.max(.02,...ratios.map(r=>1-r));
  const X=i=>48+800*i/Math.max(1,ratios.length-1),Y=r=>22+(1-r)/maxDrop*61;
  const svg=$('chart'),ns='http://www.w3.org/2000/svg';svg.replaceChildren();
  const el=(name,attrs,text)=>{const n=document.createElementNS(ns,name);for(const [k,v]of Object.entries(attrs))n.setAttribute(k,v);if(text)n.textContent=text;svg.append(n);return n;};
  el('line',{x1:48,y1:22,x2:848,y2:22,stroke:'#d4dfe1','stroke-dasharray':'4 4'});
  el('text',{x:4,y:25,fill:'#85989f','font-size':10},'0%');
  el('text',{x:4,y:85,fill:'#85989f','font-size':10},`−${(100*maxDrop).toFixed(1)}%`);
  const path=ratios.map((r,i)=>`${i?'L':'M'}${X(i)},${Y(r)}`).join(' ');
  el('path',{d:`${path} L${X(ratios.length-1)},22 L48,22 Z`,fill:accent(),'fill-opacity':.06});
  el('path',{d:path,fill:'none',stroke:accent(),'stroke-width':2});
  const maxVisible=30,stride=Math.max(1,Math.ceil(ratios.length/maxVisible));
  for(let i=0;i<ratios.length;i++){
    if(i%stride&&i!==state.index&&i!==ratios.length-1)continue;
    const chosen=i===state.index&&!state.final;
    const c=el('circle',{cx:X(i),cy:Y(ratios[i]),r:chosen?5.5:3,fill:chosen?accent():'#fff',stroke:accent(),'stroke-width':1.5,class:'milestone',role:'button',tabindex:0,'aria-label':`${i}단계, ${percent(ratios[i])}% 감소`});
    c.addEventListener('click',()=>setStep(i));c.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();setStep(i);}});
    if(chosen){el('text',{x:Math.max(80,Math.min(815,X(i))),y:Math.max(12,Y(ratios[i])-10),'text-anchor':'middle',fill:accent(),'font-size':11,'font-weight':600},`${i}단계 · −${percent(ratios[i])}%`);}
  }
  el('text',{x:48,y:103,fill:'#93a2a7','font-size':9},'기준');
  el('text',{x:848,y:103,fill:'#93a2a7','font-size':9,'text-anchor':'end'},`${series.stages.length}단계`);
}

function drawSection(){
  if(!data)return;
  const canvas=$('section'),rect=canvas.getBoundingClientRect(),dpr=Math.min(devicePixelRatio||1,2);
  canvas.width=Math.round(rect.width*dpr);canvas.height=Math.round(rect.height*dpr);
  const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);const w=rect.width,h=rect.height;
  ctx.fillStyle='#fff';ctx.fillRect(0,0,w,h);
  const bounds=state.zoom?(state.source==='IMP'?[.97,1.19,-1.38,-1.10]:[1.10,1.47,-1.44,-1.24]):[.97,1.73,-1.44,-1.08];
  const [xmin,xmax,ymin,ymax]=bounds,pad={left:48,right:15,top:25,bottom:39};
  const scale=Math.min((w-pad.left-pad.right)/(xmax-xmin),(h-pad.top-pad.bottom)/(ymax-ymin));
  const ox=pad.left+(w-pad.left-pad.right-scale*(xmax-xmin))/2,oy=pad.top+(h-pad.top-pad.bottom-scale*(ymax-ymin))/2;
  const X=x=>ox+(x-xmin)*scale,Y=y=>oy+(ymax-y)*scale;
  sectionFrame={scale,xmin,xmax,ymin,ymax,ox,oy,
    wall:segmentsForWall(shownWall(),false).map(([a,b])=>[...a,...b]),
    ghost:record&&$('display').value==='overlay'?segmentsForWall(data.baseline.wall,false).map(([a,b])=>[...a,...b]):[]};
  ctx.font='9px Segoe UI, sans-serif';ctx.lineWidth=1;
  const dx=state.zoom ? .05 : .1;
  for(let x=Math.ceil(xmin/dx)*dx;x<xmax;x+=dx){ctx.strokeStyle='#eef2f3';ctx.beginPath();ctx.moveTo(X(x),Y(ymin));ctx.lineTo(X(x),Y(ymax));ctx.stroke();ctx.fillStyle='#84989f';ctx.textAlign='center';ctx.fillText(x.toFixed(2),X(x),Y(ymin)+17);}
  for(let y=Math.ceil(ymin/.05)*.05;y<ymax;y+=.05){ctx.strokeStyle='#eef2f3';ctx.beginPath();ctx.moveTo(X(xmin),Y(y));ctx.lineTo(X(xmax),Y(y));ctx.stroke();ctx.fillStyle='#84989f';ctx.textAlign='right';ctx.fillText(y.toFixed(2),X(xmin)-7,Y(y)+3);}
  ctx.fillStyle='#697e85';ctx.textAlign='center';ctx.fillText('R [m]',(X(xmin)+X(xmax))/2,h-6);ctx.save();ctx.translate(12,h/2);ctx.rotate(-Math.PI/2);ctx.fillText('Z [m]',0,0);ctx.restore();
  ctx.save();ctx.beginPath();ctx.rect(X(xmin),Y(ymax),scale*(xmax-xmin),scale*(ymax-ymin));ctx.clip();
  function line(points,color,width,dash=[]){ctx.beginPath();ctx.strokeStyle=color;ctx.lineWidth=width;ctx.setLineDash(dash);points.forEach((p,i)=>{if(i)ctx.lineTo(X(p[0]),Y(p[1]));else ctx.moveTo(X(p[0]),Y(p[1]));});ctx.stroke();ctx.setLineDash([]);}
  const display=$('display').value;
  if(display!=='candidate'||!record)line(data.baseline.wall,'#96a5ab',2,record?[5,4]:[]);
  if(record&&display!=='baseline'){
    line(record.wall,'#626f75',2.2);
    for(const curve of Object.values(record.active_curves??{}))line(curve,$('heat').checked?'#899da7':accent(),3.4);
  }
  if(heatVisible()){
    for(const [r1,z1,r2,z2,q]of state.profile.segments){const rgb=heatRGB(q,heatMax()).map(v=>Math.round(v*255));line([[r1,z1],[r2,z2]],`rgb(${rgb.join(',')})`,4.2);}
    const [r,z]=state.profile.peak_point,x=X(r),y=Y(z);
    ctx.fillStyle='#fff';ctx.strokeStyle='#283e46';ctx.lineWidth=1.5;ctx.beginPath();ctx.moveTo(x,y-5);ctx.lineTo(x+5,y);ctx.lineTo(x,y+5);ctx.lineTo(x-5,y);ctx.closePath();ctx.fill();ctx.stroke();
  }
  for(const {legs:pair,ghost} of visibleLegs())for(const [name,leg]of Object.entries(pair)){
    if(!leg.resolved)continue;
    line(leg.points,legColors[name],ghost?1.4:2,ghost?[4,4]:[]);
    const [r,z]=leg.strike_RZ_m;
    ctx.beginPath();ctx.arc(X(r),Y(z),ghost?5:3.5,0,2*Math.PI);
    ctx.fillStyle=ghost?'#fff':legColors[name];ctx.strokeStyle=legColors[name];ctx.lineWidth=1.5;ctx.fill();ctx.stroke();
  }
  if(visibleLegs().length){
    const [r,z]=legs.xpoint_RZ_m,x=X(r),y=Y(z);ctx.strokeStyle='#253e47';ctx.lineWidth=1.7;
    ctx.beginPath();ctx.moveTo(x-4,y-4);ctx.lineTo(x+4,y+4);ctx.moveTo(x+4,y-4);ctx.lineTo(x-4,y+4);ctx.stroke();
    ctx.fillStyle='#253e47';ctx.textAlign='left';ctx.font='10px Segoe UI, sans-serif';ctx.fillText('X-point',x+8,y-6);
  }
  ctx.restore();
}

function initScene(){
  try{
    renderer=new THREE.WebGLRenderer({antialias:true,alpha:true,preserveDrawingBuffer:true});
    renderer.setPixelRatio(Math.min(devicePixelRatio||1,1.75));renderer.outputColorSpace=THREE.SRGBColorSpace;
    renderer.setClearColor(0xffffff,0);$('scene').appendChild(renderer.domElement);
    scene=new THREE.Scene();camera=new THREE.PerspectiveCamera(37,1,.01,40);
    scene.add(new THREE.HemisphereLight(0xffffff,0xa1b9c1,2.5));
    const light=new THREE.DirectionalLight(0xffffff,2.7);light.position.set(3,6,4);scene.add(light);
    controls=new OrbitControls(camera,renderer.domElement);controls.enableDamping=false;controls.minDistance=.2;controls.maxDistance=12;controls.addEventListener('change',()=>{hideProbe();render();});
    group=new THREE.Group();scene.add(group);preset('focus');resizeScene();
  }catch(e){renderer=null;$('webgl-fallback').hidden=false;$('export').disabled=true;console.warn('WebGL unavailable:',e.message);}
}
function clearGroup(){
  pickWall=null;
  if(!group)return;
  group.traverse(o=>{o.geometry?.dispose();if(o.material)for(const m of Array.isArray(o.material)?o.material:[o.material])m.dispose();});group.clear();
}
function buildScene(){
  hideProbe();if(!renderer||!data)return;clearGroup();
  const lower=!$('context').checked,display=$('display').value,angle=Number($('angle').value)*Math.PI/180,start=Math.PI/2+(2*Math.PI-angle)/2;
  const sectors=Math.max(12,Math.ceil(angle/(2*Math.PI)*80));
  function surface(segments,color,opacity=1,heat=false,pickable=false){
    const positions=[],colors=[],uv=[];
    for(const s of segments){
      const a=s[0],b=s[1],rgb=heat?new THREE.Color().setRGB(...heatRGB(s[2],heatMax()),THREE.SRGBColorSpace).toArray():null;
      for(let j=0;j<sectors;j++){
        const p0=revolve(a,start+angle*j/sectors),p1=revolve(b,start+angle*j/sectors),p2=revolve(b,start+angle*(j+1)/sectors),p3=revolve(a,start+angle*(j+1)/sectors);
        positions.push(...p0,...p1,...p2,...p0,...p2,...p3);
        if(pickable)uv.push(...segmentTriangleRZ(a,b));
        if(heat)for(let k=0;k<6;k++)colors.push(...rgb);
      }
    }
    if(!positions.length)return;
    const geom=new THREE.BufferGeometry();geom.setAttribute('position',new THREE.Float32BufferAttribute(positions,3));
    if(pickable)geom.setAttribute('uv',new THREE.Float32BufferAttribute(uv,2));
    if(heat)geom.setAttribute('color',new THREE.Float32BufferAttribute(colors,3));
    geom.computeVertexNormals();
    const mat=heat?new THREE.MeshBasicMaterial({vertexColors:true,side:THREE.DoubleSide,polygonOffset:true,polygonOffsetFactor:-3,polygonOffsetUnits:-3}):new THREE.MeshStandardMaterial({color,side:THREE.DoubleSide,roughness:.72,metalness:.12,transparent:opacity<1,opacity,depthWrite:opacity===1});
    const mesh=new THREE.Mesh(geom,mat);group.add(mesh);if(pickable)pickWall=mesh;
  }
  function lines(wall,color,ghost=false){
    for(const phi of [start,start+angle,...(ghost?[]:[start+angle/2])]){
      const pos=[];for(const [a,b]of segmentsForWall(wall,lower))pos.push(...revolve(a,phi),...revolve(b,phi));
      const geom=new THREE.BufferGeometry();geom.setAttribute('position',new THREE.Float32BufferAttribute(pos,3));
      const mat=ghost?new THREE.LineDashedMaterial({color,transparent:true,opacity:.7,dashSize:.013,gapSize:.009,depthTest:false}):new THREE.LineBasicMaterial({color,transparent:true,opacity:.8});
      const line=new THREE.LineSegments(geom,mat);line.computeLineDistances();line.renderOrder=ghost?4:2;group.add(line);
    }
  }
  const shown=shownWall();
  surface(segmentsForWall(shown,lower),'#b9c8cc',1,false,true);
  if(record&&display!=='baseline')for(const curve of Object.values(record.active_curves??{}))surface(segmentsForWall(curve,lower),$('heat').checked?'#a7bac3':accent());
  lines(shown,'#425f68');
  if(record&&display==='overlay')lines(data.baseline.wall,'#6e8090',true);
  if($('plasma').checked){surface(segmentsForWall(data.baseline.lcfs,false),'#d4aac0',.13);lines(data.baseline.lcfs,'#bf8ba6');}
  if(heatVisible()){
    // Pick the colored receiving cells themselves, including cut-edge cells.
    surface(state.profile.segments.filter(([r,z,r2,z2])=>!lower||(z<-1.05&&z2<-1.05)).map(([r,z,r2,z2,q])=>[[r,z],[r2,z2],q]),null,1,true,true);
    const p=revolve(state.profile.peak_point,start+.02),mesh=new THREE.Mesh(new THREE.SphereGeometry(.007,12,8),new THREE.MeshBasicMaterial({color:'#ffffff',depthTest:false}));mesh.position.set(...p);mesh.renderOrder=5;group.add(mesh);
  }
  const layers=visibleLegs();
  const marker=(point,phi,color,radius)=>{
    const m=new THREE.Mesh(new THREE.SphereGeometry(radius,12,8),new THREE.MeshBasicMaterial({color,depthTest:false}));
    m.position.set(...revolve(point,phi));m.renderOrder=8;group.add(m);
  };
  for(const {legs:pair,ghost}of layers)for(const [name,leg]of Object.entries(pair)){
    if(!leg.resolved)continue;
    // The translucent swept surface is axisymmetric; the colored curves are poloidal cuts.
    if(!ghost)surface(segmentsForWall(leg.points,false),legColors[name],.18);
    for(const phi of [start,start+angle]){
      const geometry=new THREE.BufferGeometry().setFromPoints(leg.points.map(p=>new THREE.Vector3(...revolve(p,phi))));
      const material=ghost?new THREE.LineDashedMaterial({color:legColors[name],dashSize:.009,gapSize:.006,transparent:true,opacity:.65,depthTest:false}):new THREE.LineBasicMaterial({color:legColors[name],depthTest:false});
      const curve=new THREE.Line(geometry,material);curve.computeLineDistances();curve.renderOrder=ghost?6:7;group.add(curve);
      if(ghost){marker(leg.strike_RZ_m,phi,legColors[name],.0048);marker(leg.strike_RZ_m,phi,'#ffffff',.003);}
      else marker(leg.strike_RZ_m,phi,legColors[name],.0035);
    }
  }
  if(layers.length)for(const phi of [start,start+angle]){
    const [r,z]=legs.xpoint_RZ_m,points=[[r-.004,z-.004],[r+.004,z+.004],[r-.004,z+.004],[r+.004,z-.004]];
    const geometry=new THREE.BufferGeometry().setFromPoints(points.map(p=>new THREE.Vector3(...revolve(p,phi))));
    const cross=new THREE.LineSegments(geometry,new THREE.LineBasicMaterial({color:'#253e47',depthTest:false}));cross.renderOrder=9;group.add(cross);
  }
  render();
}
function render(){if(renderer&&camera)renderer.render(scene,camera);}
function resizeScene(){hideProbe();if(!renderer)return;const r=$('scene').getBoundingClientRect();renderer.setSize(r.width,r.height);camera.aspect=r.width/r.height;camera.updateProjectionMatrix();render();}
function preset(kind){
  if(!renderer)return;
  if(kind==='overview'){camera.position.set(4.2,2.8,5.8);controls.target.set(0,-.1,0);}
  else if(kind==='top'){camera.position.set(.05,4.8,.1);controls.target.set(0,-1.25,0);}
  else {camera.position.set(2.9,.75,4.15);controls.target.set(0,-1.25,0);}
  controls.update();render();
}
function setStep(i){stop();state.final=false;state.index=Math.min(series.stages.length,Math.max(0,i));update().catch(error);}

function hideProbe(){
  probePointer=null;cancelAnimationFrame(probeFrame);probeFrame=0;$('q-tooltip').hidden=true;
}
function showProbe(point,event,ghost=false){
  const hit=!ghost&&state.profile?nearestHeatSegment(point,state.profile.segments,2e-5):null;
  const name=ghost||!record||$('display').value==='baseline'?'기준 형상':'선택 형상';
  $('probe-title').textContent=`${name} · ${state.source}`;
  $('probe-value').textContent=hit?`q = ${formatHeatQ(hit.q)} MW/m²`:'q · 기록 없음';
  $('probe-position').textContent=`R ${point[0].toFixed(5)} m · Z ${point[1].toFixed(5)} m`;
  $('probe-detail').textContent=ghost?'기준 형상만 선택하면 기준 q를 확인할 수 있습니다.':hit?`${state.profile.donors??1280} 경로 · 표시 셀의 면적 가중 평균`:profileLoading?'해당 형상의 저장 분포를 불러오는 중입니다.':state.profile?'이 위치는 열분포 격자 범위 밖입니다.':'해당 형상의 열분포를 불러오지 못했습니다.';
  $('probe-swatch').hidden=!hit;
  if(hit)$('probe-swatch').style.background=`rgb(${heatRGB(hit.q,heatMax()).map(v=>Math.round(v*255)).join(',')})`;
  const box=$('q-tooltip');box.hidden=false;
  box.style.left=Math.max(8,Math.min(event.clientX+15,innerWidth-box.offsetWidth-8))+'px';
  box.style.top=Math.max(8,Math.min(event.clientY+15,innerHeight-box.offsetHeight-8))+'px';
}
function probeSection(event){
  if(!sectionFrame)return;
  const f=sectionFrame,rect=$('section').getBoundingClientRect();
  const point=[f.xmin+(event.clientX-rect.left-f.ox)/f.scale,f.ymax-(event.clientY-rect.top-f.oy)/f.scale];
  if(point[0]<f.xmin||point[0]>f.xmax||point[1]<f.ymin||point[1]>f.ymax){hideProbe();return;}
  const main=nearestHeatSegment(point,f.wall,7/f.scale),ghost=nearestHeatSegment(point,f.ghost,7/f.scale);
  if(ghost&&(!main||ghost.distance+1/f.scale<main.distance))showProbe(ghost.point,event,true);
  else if(main)showProbe(main.point,event);else hideProbe();
}
function probeScene(event){
  if(!renderer||!pickWall)return;
  const rect=renderer.domElement.getBoundingClientRect();
  probeRay.setFromCamera(new THREE.Vector2(2*(event.clientX-rect.left)/rect.width-1,1-2*(event.clientY-rect.top)/rect.height),camera);
  const hit=probeRay.intersectObject(pickWall,false)[0];
  if(hit?.uv)showProbe([hit.uv.x,hit.uv.y],event);else hideProbe();
}
function bindProbe(element,probe){
  element.addEventListener('pointermove',e=>{
    if(e.buttons){hideProbe();return;}
    probePointer={clientX:e.clientX,clientY:e.clientY};
    if(!probeFrame)probeFrame=requestAnimationFrame(()=>{probeFrame=0;const p=probePointer;if(p)probe(p);});
  });
  let down=null;
  element.addEventListener('pointerdown',e=>{down=[e.clientX,e.clientY];hideProbe();});
  element.addEventListener('pointerup',e=>{if(down&&Math.hypot(e.clientX-down[0],e.clientY-down[1])<4)probe(e);down=null;});
  element.addEventListener('pointerleave',hideProbe);
}

function bind(){
  bindProbe($('section'),probeSection);bindProbe($('scene'),probeScene);
  addEventListener('scroll',hideProbe,{passive:true});
  document.addEventListener('keydown',e=>{if(e.key==='Escape')hideProbe();});
  $('sources').addEventListener('click',e=>{const b=e.target.closest('button');if(b){state.source=b.dataset.source;selectSeries();}});
  $('radii').addEventListener('click',e=>{const b=e.target.closest('button');if(b){state.radius=Number(b.dataset.radius);selectSeries();}});
  $('runs').addEventListener('change',()=>selectSeries($('runs').value));
  $('timeline').addEventListener('input',()=>setStep(Number($('timeline').value)));
  $('stage-select').addEventListener('change',()=>setStep(Number($('stage-select').value)));
  $('prev').addEventListener('click',()=>setStep(state.final?(record?.archived_stage??series.stages.length):state.index-1));
  $('next').addEventListener('click',()=>setStep(state.final&&!record?.additional_validation?series.stages.length:state.index+1));
  $('final').addEventListener('click',()=>{stop();state.final=true;if(series.final?.additional_validation)state.index=series.final.archived_stage;update().catch(error);});
  $('play').addEventListener('click',()=>{
    if(timer){stop();return;}
    if(state.final||state.index===series.stages.length){state.index=0;state.final=false;update().catch(error);}
    $('play').textContent='Ⅱ';$('play').setAttribute('aria-label','자동 재생 일시정지');
    timer=setInterval(()=>{if(state.index>=series.stages.length){stop();return;}state.index++;update().catch(error);},900);
  });
  for(const id of ['display','context','plasma','heat','legs'])$(id).addEventListener('change',()=>update().catch(error));
  $('angle').addEventListener('input',()=>{$('angle-value').textContent=$('angle').value+'°';buildScene();});
  $('section-zoom').addEventListener('click',()=>{state.zoom=!state.zoom;update().catch(error);});
  for(const id of ['focus','overview','top'])$(id).addEventListener('click',()=>preset(id));
  $('export').addEventListener('click',()=>{
    stop();render();const a=$('image-download');
    a.download=`D3D-${state.source}-R${state.radius}-${state.final?'refined':'step'+state.index}.png`;
    a.href=renderer.domElement.toDataURL('image/png');$('image-preview').src=a.href;$('image-dialog').showModal();
  });
  new ResizeObserver(()=>{drawSection();resizeScene();}).observe($('section'));
  document.addEventListener('visibilitychange',()=>{if(document.hidden)stop();});
}

try{
  const response=await fetch('data.json');if(!response.ok)throw new Error('data.json을 읽을 수 없습니다. launch.ps1로 실행해 주세요.');
  const body=await response.text();data=JSON.parse(body);
  const dataHash=await digestText(body);
  const [stageResponse,savedResponse]=await Promise.all([fetch('stage-profiles.json'),fetch('full-saved-profiles.json')]);
  if(!stageResponse.ok||!savedResponse.ok)throw new Error('단계별 열분포 색인을 읽을 수 없습니다.');
  stages=validateStageIndex(await stageResponse.json(),data,dataHash);
  savedProfiles=validateSavedIndex(await savedResponse.json(),data,dataHash);
  const qualifiedResponse=await fetch('qualified-finals.json');
  if(qualifiedResponse.ok){
    const index=await qualifiedResponse.json(),evidence={};
    if(!Array.isArray(index.entries))throw new Error('추가 정밀 결과 색인 형식 오류');
    await Promise.all(index.entries.map(async entry=>{
      if(!/^qualified-evidence\/[a-f0-9]{64}\.json$/.test(entry.evidence?.path??'')||!/^[a-f0-9]{64}$/.test(entry.evidence?.sha256??''))throw new Error('추가 검증 파일 경로/해시 오류');
      const r=await fetch(entry.evidence.path);if(!r.ok)throw new Error('추가 검증 근거 파일 누락');
      const text=await r.text();if(await digestText(text)!==entry.evidence.sha256)throw new Error('추가 검증 근거 해시 불일치');
      evidence[entry.final.id]=JSON.parse(text);
    }));
    const extended=applyQualifiedFinals(index,data,savedProfiles,dataHash,evidence);
    data=extended.data;savedProfiles=extended.savedProfiles;
  }else if(qualifiedResponse.status!==404)throw new Error('추가 정밀 결과 읽기 실패');
  const distributionResponse=await fetch('distribution-finals.json');
  if(distributionResponse.ok){
    const index=await distributionResponse.json(),evidence={};
    if(index.policy!==distributionPolicy||!Array.isArray(index.entries)||!Array.isArray(index.reassessment))throw new Error('분포 수렴 검증 색인 오류');
    for(const entry of index.entries){
      if(!/^distribution-evidence\/[a-f0-9]{64}\.json$/.test(entry.evidence?.path??''))throw new Error('분포 검증 경로 오류');
      const response=await fetch(entry.evidence.path);if(!response.ok)throw new Error('분포 검증 근거 누락');
      const text=await response.text();if(await digestText(text)!==entry.evidence.sha256)throw new Error('분포 검증 해시 불일치');
      evidence[entry.final.id]=JSON.parse(text);
    }
    const extended=applyQualifiedFinals(index,data,savedProfiles,dataHash,evidence);
    data=applyReassessment(extended.data,index.reassessment);savedProfiles=extended.savedProfiles;
  }else if(distributionResponse.status!==404)throw new Error('분포 검증 읽기 실패');
  try{
    const res=await fetch('legs.json');if(!res.ok)throw new Error('leg 파일 읽기 실패');
    const loaded=await res.json();
    const digest=await crypto.subtle.digest('SHA-256',new TextEncoder().encode(body));
    const hash=Array.from(new Uint8Array(digest),v=>v.toString(16).padStart(2,'0')).join('');
    if(loaded.schema!=='d3d-viewer-poloidal-legs.v1'||loaded.data_sha256!==hash)throw new Error('형상과 leg 파일 버전 불일치');
    legs=loaded;
  }catch(e){legError='Leg를 표시할 수 없습니다: '+e.message;$('legs').checked=false;$('legs').disabled=true;console.warn(legError);}
  initScene();bind();selectSeries();
}catch(e){error(e);}
