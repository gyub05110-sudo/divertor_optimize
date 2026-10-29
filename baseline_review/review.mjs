import * as THREE from 'three';
import {OrbitControls} from '../input/D3D/optimized/d3d-divertor-floor-shelf/viewer3d/vendor/OrbitControls.js';
const d=await fetch('baseline_review.json').then(r=>{if(!r.ok)throw Error('기준 형상 자료 로드 실패');return r.json()});
const $=id=>document.getElementById(id),host=$('scene');
const scene=new THREE.Scene();scene.background=new THREE.Color('#f8fafc');
const camera=new THREE.PerspectiveCamera(35,1,.01,50);
const renderer=new THREE.WebGLRenderer({antialias:true});renderer.setPixelRatio(Math.min(devicePixelRatio,1.7));host.append(renderer.domElement);
const controls=new OrbitControls(camera,renderer.domElement);controls.enableDamping=false;
let group,resultData,sectionMode=true;const colors={IMP:'#008d99',OMP:'#36a25b',inner:'#bc35a2',outer:'#d34730'};
function point(p,phi=0){return new THREE.Vector3(p[0]*Math.cos(phi),p[1],p[0]*Math.sin(phi));}
function line(points,color,phi=0,opacity=1){if(sectionMode&&phi!==0)return;const g=new THREE.BufferGeometry().setFromPoints(points.map(p=>point(p,phi)));group.add(new THREE.Line(g,new THREE.LineBasicMaterial({color,transparent:opacity<1,opacity,depthTest:false})));}
function surface(points,color,opacity,sectors=52){
 if(sectionMode)return;
 const end=+$('angle').value*Math.PI/180,vertices=[];
 for(let i=1;i<points.length;i++)for(let j=0;j<sectors;j++){
  const a=point(points[i-1],end*j/sectors),b=point(points[i],end*j/sectors),c=point(points[i],end*(j+1)/sectors),e=point(points[i-1],end*(j+1)/sectors);
  for(const v of [a,b,c,a,c,e])vertices.push(v.x,v.y,v.z);
 }
 const g=new THREE.BufferGeometry();g.setAttribute('position',new THREE.Float32BufferAttribute(vertices,3));g.computeVertexNormals();
 const mesh=new THREE.Mesh(g,new THREE.MeshBasicMaterial({color,opacity,transparent:opacity<1,side:THREE.DoubleSide,depthWrite:opacity>=1}));group.add(mesh);
}
function dot(p,color,r=.006){const mesh=new THREE.Mesh(new THREE.SphereGeometry(r,14,10),new THREE.MeshBasicMaterial({color,depthTest:false}));mesh.position.copy(point(p));group.add(mesh);}
function label(p,text,color){
 const c=document.createElement('canvas');let ctx=c.getContext('2d');ctx.font='bold 27px Malgun Gothic, sans-serif';c.width=Math.ceil(ctx.measureText(text).width)+24;c.height=54;
 ctx=c.getContext('2d');ctx.font='bold 27px Malgun Gothic, sans-serif';ctx.fillStyle='#fffffff0';ctx.fillRect(0,0,c.width,c.height);ctx.fillStyle=color;ctx.fillText(text,12,37);
 const m=new THREE.SpriteMaterial({map:new THREE.CanvasTexture(c),depthTest:false});const s=new THREE.Sprite(m);s.position.copy(point(p));s.scale.set(.045*c.width/54,.045,1);group.add(s);
}
function pfr(phi){const poly=d.pfr_polygon,s=new THREE.Shape();s.moveTo(...poly[0]);for(const p of poly.slice(1))s.lineTo(...p);s.closePath();
 if(sectionMode&&phi!==0)return;
 const g=new THREE.ShapeGeometry(s),positions=g.attributes.position;
 for(let i=0;i<positions.count;i++){const r=positions.getX(i),z=positions.getY(i);positions.setXYZ(i,r*Math.cos(phi),z,r*Math.sin(phi));}
 group.add(new THREE.Mesh(g,new THREE.MeshBasicMaterial({color:'#dba4c5',opacity:.55,transparent:true,side:THREE.DoubleSide,depthWrite:false})));
}
function rebuild(){
 if(group){scene.remove(group);group.traverse(o=>{o.geometry?.dispose();if(o.material){o.material.map?.dispose();o.material.dispose()}})}group=new THREE.Group();scene.add(group);
 const end=+$('angle').value*Math.PI/180;
 if($('wall').checked){surface(d.wall,'#9caab7',.17);line(d.wall,'#7d8fa0');line(d.wall,'#9caab7',end,.5);}
 if($('lcfs').checked){surface(d.lcfs,'#4f77a5',.06);line(d.lcfs,'#4f77a5',0,.7);}
 if($('design').checked)for(const seg of d.existing_design_segments){surface(d.wall.slice(seg,seg+2),'#ed9c18',.50);line(d.wall.slice(seg,seg+2),'#ed9c18');}
 if($('pfr').checked){pfr(0);if(end<Math.PI*2)pfr(end);label([1.18,-1.385],'PFR · 열분포 미정','#95547f');}
 for(const src of ['IMP','OMP'])if($(src.toLowerCase()).checked){
  for(const r of d.sources[src].paths)line(r.points,colors[src],0,.55);
  const regions=resultData?.boundary_mode==='continuous'&&resultData.design_surfaces
    ?resultData.design_surfaces.filter(s=>(src==='IMP'?[69,70]:[58,60]).includes(s.original_segment)).map(s=>({points:s.endpoints_RZ_m}))
    :d.optimization_footprints.filter(r=>r.source===src);
  for(const r of regions){surface(r.points,colors[src],.6);line(r.points,colors[src]);}
  dot(d.sources[src].midplane_RZ_m,colors[src],.009);label([d.sources[src].midplane_RZ_m[0],d.sources[src].midplane_RZ_m[1]+.04],src+' 열원',colors[src]);
 }
 if($('imp').checked){dot([1.01748898,-1.15],'#008d99',.004);label([1.10,-1.13],'IMP: Z=-1.15부터 아래','#008d99');}
 if($('legs').checked){for(const name of ['inner','outer']){const p=d.legs[name].points;surface(p,colors[name],.18);line(p,colors[name]);line(p,colors[name],end,.8);dot(d.legs[name].strike_RZ_m,colors[name]);
  const curve=new THREE.CatmullRomCurve3(p.map(x=>point(x)));group.add(new THREE.Mesh(new THREE.TubeGeometry(curve,p.length,.0013,6,false),new THREE.MeshBasicMaterial({color:colors[name],depthTest:false})));}
  dot(d.xpoint_RZ_m,'#1d2e40',.008);label([1.17659,-1.215],'X-point','#1d2e40');label([1.09,-1.257],'inner leg',colors.inner);label([1.30,-1.32],'outer leg',colors.outer);
 }
 if(resultData&&$('resultShow').checked){
  surface(resultData.wall,'#6956b3',.28);line(resultData.wall,'#6956b3');
  const h=resultData.heat[$('heatWidth').value],peak=Math.max(...h.q_over_P);
  h.points.forEach((p,i)=>dot(p,new THREE.Color().setHSL(.66*(1-h.q_over_P[i]/(peak||1)),1,.48),.0025));
  const audit=resultData.validity_audit_summary;
  const validation=audit?.verification_status==='passed'?'형상·해상도 검증 통과':audit?.verification_status?'추가 확인 필요':resultData.numerically_verified_feasible?'형상 제약 재검증 통과':'수치 재검증 전';
  $('resultStatus').textContent=`${resultData.source} · ${resultData.replay_only?'기존 후보 · 수정 모델 재평가':resultData.smoke_test?'소규모 연결 시험':resultData.qualified_search?'본 탐색':'기준 계산'} · ${validation} · 최대 q/P=${peak.toFixed(4)} m⁻². 파랑→빨강: 0→최대. PFR 열분포는 미정. ${resultData.validity_note||''}`;
 }
 $('angleValue').textContent=$('angle').value+'°';render();
}
function render(){renderer.render(scene,camera)}
function view(name){sectionMode=name==='section';camera.up.set(0,1,0);if(name==='section'){camera.position.set(1.41,-1.09,1.6);controls.target.set(1.41,-1.09,0)}
 else if(name==='focus'){camera.position.set(2.8,-.2,1.8);controls.target.set(1.16,-1.24,0)}
 else{camera.position.set(5,2.3,5.5);controls.target.set(0,-.32,0)}controls.update();rebuild();render();}
function resize(){const w=host.clientWidth,h=host.clientHeight;renderer.setSize(w,h,false);camera.aspect=w/h;camera.updateProjectionMatrix();render();}
controls.addEventListener('change',render);new ResizeObserver(resize).observe(host);
for(const id of ['wall','lcfs','legs','design','imp','omp','pfr'])$(id).addEventListener('change',rebuild);
$('angle').addEventListener('input',rebuild);for(const id of ['section','focus','overview'])$(id).onclick=()=>view(id);
$('resultShow').addEventListener('change',rebuild);$('heatWidth').addEventListener('change',rebuild);
$('resultFile').addEventListener('change',async()=>{
 try{
  const value=JSON.parse(await $('resultFile').files[0].text());
  const pointsValid=ps=>Array.isArray(ps)&&ps.length>0&&ps.length<100000&&ps.every(p=>Array.isArray(p)&&p.length===2&&p.every(Number.isFinite));
  if(value.schema!=='independent-midplane-viewer.v1'||!['IMP','OMP'].includes(value.source)||!pointsValid(value.wall))throw Error('독립 최적화 결과 형식이 아닙니다');
  for(const width of ['1','2','3','5','10']){const h=value.heat?.[width];if(!h||!pointsValid(h.points)||h.points.length!==h.q_over_P.length||!h.q_over_P.every(q=>Number.isFinite(q)&&q>=0))throw Error('열유속 데이터 오류');}
  resultData=value;rebuild();
 }catch(err){$('resultStatus').textContent=err.message;}
});
const rows=[['X-point',d.xpoint_RZ_m],['inner leg → 벽',d.legs.inner.strike_RZ_m],['outer leg → 벽',d.legs.outer.strike_RZ_m],['IMP separatrix / 중간면',d.sources.IMP.midplane_RZ_m],['OMP separatrix / 중간면',d.sources.OMP.midplane_RZ_m]];
for(const [name,p]of rows){const tr=document.createElement('tr');for(const text of [name,`R=${p[0].toFixed(6)}, Z=${p[1].toFixed(6)} m`]){const td=document.createElement('td');td.textContent=text;tr.append(td)}$('positions').append(tr)}
for(const src of ['IMP','OMP']){const tr=document.createElement('tr');tr.innerHTML=`<td>${src} 최초 충돌 원본 segment</td><td>${Object.keys(d.sources[src].hit_counts).join(', ')}</td>`;$('positions').append(tr)}
for(const [name,value]of [['IMP 설계영역','Z<=-1.15 m: segment 70 하부 + 69 일부'],['OMP 제외 영역','바깥쪽 경사 벽 segment 53 / 54 / 55: 고정, 열부하 계산 유지']]){const tr=document.createElement('tr');for(const text of [name,value]){const td=document.createElement('td');td.textContent=text;tr.append(td)}$('positions').append(tr)}
rebuild();resize();view('section');$('status').textContent='기준 형상 로드 완료 · 실제 IMP/OMP의 하부 최초 충돌 · 각 800개 표본 · 두 leg의 ψX 등고선 확인';
