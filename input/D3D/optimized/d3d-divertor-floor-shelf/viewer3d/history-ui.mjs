import {casesForFile,fileSHA256} from './registry.mjs';
import {deviceOf,casesForDevice,chooseInitialCase,stageDescription} from './device.mjs';
const $=id=>document.getElementById(id);
async function read(url){const r=await fetch('./'+url);if(!r.ok)throw Error('결과 데이터 로드 실패');return r.json();}

export function setupHistory({catalog,initialEntry,initialBundle,initialStage,applyStage,unavailableView}){
 let entry=initialEntry,bundle=initialBundle,index=bundle.stages.length-1,cache,epoch=0,available=true,device=deviceOf(entry);
 const remembered=new Map([[device,entry.id]]);
 function seed(){cache=new Map([[bundle.stages[0].id,bundle.models.baseline],[bundle.stages.at(-1).id,bundle.models.candidate]]);}
 seed();
 function fillCases(entries){
  $('equilibrium').replaceChildren();
  for(const c of entries){const option=document.createElement('option');option.value=c.id;
   option.textContent=deviceOf(c)==='cpd2'?`${c.run_name} · 고정 평형`:`${c.shot_name} · ${c.time_ms} ms · ${c.run_name}`;$('equilibrium').append(option);}
  if(entries.some(c=>c.id===entry.id))$('equilibrium').value=entry.id;
 }
 function fillStages(){
  $('stageSelect').replaceChildren();
  for(const s of bundle.stages){const option=document.createElement('option');option.value=s.id;option.textContent=s.label+(s.refinement_rejected?' [재검증 탈락]':'');$('stageSelect').append(option);}
  $('stageSelect').disabled=false;$('stageRange').disabled=false;$('stageRange').max=bundle.stages.length-1;
 }
 function missing(message){
  available=false;epoch++;unavailableView(message);$('fileMessage').textContent=message;
  $('stageSelect').replaceChildren();$('stageSelect').disabled=true;$('stageRange').disabled=true;
  $('previousStage').disabled=true;$('nextStage').disabled=true;$('stageInfo').textContent='연결된 최적화 단계 없음';
 }
 async function loadStage(nextIndex,viewParams){
  if(!available||nextIndex<0||nextIndex>=bundle.stages.length)return;
  const token=++epoch,s=bundle.stages[nextIndex];
  $('status').textContent='선택 단계 형상을 불러오는 중…';$('status').hidden=false;
  try{
   let detail=cache.get(s.id);
   if(!detail){detail=await read(s.detail_url);if(token!==epoch)return;cache.set(s.id,detail);}
   if(token!==epoch)return;index=nextIndex;applyStage(bundle,detail,entry,s,viewParams);
   $('stageSelect').value=s.id;$('stageRange').value=index;
   $('previousStage').disabled=index===0;$('nextStage').disabled=index===bundle.stages.length-1;
   $('stageInfo').textContent=stageDescription(bundle,s,index);
   $('stageInfo').classList.toggle('warning',s.refinement_rejected===true);
  }catch(error){if(token===epoch)missing('단계 데이터를 불러오지 못했습니다. 등록 목록에서 다시 선택해 주세요.');console.error(error);}
 }
 async function loadCase(next,stageId,viewParams){
  const token=++epoch;available=false;unavailableView('선택 평형과 이력을 불러오는 중…');
  device=deviceOf(next);$('device').value=device;
  try{
   const result=next.id===entry.id?bundle:await read(next.bundle_url);
   if(token!==epoch)return;
   entry=next;bundle=result;remembered.set(device,entry.id);seed();available=true;fillStages();
   if(!$('equilibrium').querySelector(`option[value="${entry.id}"]`))fillCases(casesForDevice(catalog.cases,device));
   $('equilibrium').value=entry.id;
   $('fileMessage').textContent=bundle.source_model==='unified_radial_exponential_v1'?`연결 완료 · ${entry.run_name} · 검증된 기준/후보 ${bundle.stages.length}단계`:device==='cpd2'?`연결 완료 · ${entry.run_name} · 정밀 표본 ${bundle.stages.length}단계`:
    `연결 완료 · ${entry.time_ms} ms · 기준 + 개선 후보 + 최종 검증 ${bundle.stages.length}단계`;
   const stage=bundle.stages.findIndex(s=>s.id===stageId);await loadStage(stage>=0?stage:bundle.stages.length-1,viewParams);
  }catch(error){if(token===epoch)missing('이 평형의 결과 데이터를 불러오지 못했습니다.');console.error(error);}
 }
 $('device').replaceChildren();
 for(const id of [...new Set(catalog.cases.map(deviceOf))]){const option=document.createElement('option');option.value=id;option.textContent=id==='cpd2'?'CPD-2':'D3D (DIII-D)';$('device').append(option);}
 $('device').value=device;
 $('device').addEventListener('change',()=>{
  device=$('device').value;const entries=casesForDevice(catalog.cases,device);fillCases(entries);
  const next=entries.find(c=>c.id===remembered.get(device))||entries[0];if(next)loadCase(next);else missing('등록된 장치 결과가 없습니다.');
 });
 $('equilibrium').addEventListener('change',()=>{const c=casesForDevice(catalog.cases,device).find(c=>c.id===$('equilibrium').value);if(c)loadCase(c);});
 $('stageSelect').addEventListener('change',()=>loadStage(bundle.stages.findIndex(s=>s.id===$('stageSelect').value)));
 $('stageRange').addEventListener('change',()=>loadStage(+$('stageRange').value));
 $('previousStage').addEventListener('click',()=>loadStage(index-1));$('nextStage').addEventListener('click',()=>loadStage(index+1));
 $('allCases').addEventListener('click',()=>{const entries=casesForDevice(catalog.cases,device);fillCases(entries);if(entries.length)loadCase(entries.find(c=>c.id===entry.id)||entries[0]);});
 $('equilibriumFile').addEventListener('change',async()=>{
  const file=$('equilibriumFile').files[0];if(!file)return;
  missing('선택한 평형 파일과 연결된 결과를 찾는 중…');const token=epoch;
  try{
   const hash=await fileSHA256(file);if(token!==epoch)return;
   const entries=casesForFile(catalog.cases,hash);
   if(!entries.length){fillCases([]);missing(`등록된 최적화 결과가 없습니다: ${file.name}. 결과를 등록한 뒤 선택해 주세요.`);return;}
   device=deviceOf(entries[0]);$('device').value=device;fillCases(entries);
   await loadCase(entries.find(c=>c.id===entry.id)||entries[0]);
  }catch(error){if(token===epoch)missing('평형 파일을 확인하지 못했습니다. 등록 목록에서 선택해 주세요.');console.error(error);}
 });
 const restore=()=>{const viewParams=new URLSearchParams(location.hash.slice(1)),next=chooseInitialCase(catalog.cases,viewParams);
  fillCases(casesForDevice(catalog.cases,deviceOf(next)));loadCase(next,viewParams.get('stage'),viewParams);};
 window.addEventListener('hashchange',restore);window.addEventListener('popstate',restore);
 fillCases(casesForDevice(catalog.cases,device));fillStages();
 return {ready:loadCase(initialEntry,initialStage),selectModel:kind=>loadStage(kind==='baseline'?0:bundle.stages.length-1)};
}
