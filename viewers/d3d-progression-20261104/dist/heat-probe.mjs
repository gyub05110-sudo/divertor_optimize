// Pure geometry and data-selection helpers for the surface q inspector.
export function nearestHeatSegment([r,z],segments,maxDistance=Infinity){
  let best=null,bestSquared=maxDistance*maxDistance;
  for(let i=0;i<segments.length;i++){
    const [ar,az,br,bz,q]=segments[i],dr=br-ar,dz=bz-az,den=dr*dr+dz*dz;
    const t=den>0?Math.max(0,Math.min(1,((r-ar)*dr+(z-az)*dz)/den)):0;
    const point=[ar+t*dr,az+t*dz],squared=(r-point[0])**2+(z-point[1])**2;
    if(squared<=bestSquared&&(!best||squared<bestSquared)){
      bestSquared=squared;best={index:i,t,point,q:q??null,distance:Math.sqrt(squared)};
    }
  }
  return best;
}
export function profileForView(baseline,source,record,display,stages=null){
  if(display==='baseline'||!record)return baseline.profiles[source]??null;
  if(record.profile)return record.profile;
  const stage=stages?.cases?.[record.id];
  return stage?.source===source?stage.profile:null;
}
export function segmentTriangleRZ(a,b){return [...a,...b,...b,...a,...b,...a];}
export function formatHeatQ(q){
  if(typeof q!=='number'||!Number.isFinite(q))return '기록 없음';
  return q>0&&q<1e-4?q.toExponential(3):q.toFixed(4);
}
export function validateStageIndex(index,data,hash){
  if(index?.schema!=='d3d-stage-profiles.v1'||!index.complete||index.data_sha256!==hash)throw new Error('단계별 열분포와 형상 파일 버전 불일치');
  for(const series of data.series)for(const row of series.stages){
    const entry=index.cases?.[row.id];
    if(!entry||entry.source!==series.source||entry.donors!==320||entry.profile!==`stage-profiles/${row.id}.json`||!validSHA(entry.sha256))throw new Error('단계별 열분포 누락 또는 후보 불일치');
  }
  return index;
}
const validSHA=value=>typeof value==='string'&&/^[a-f0-9]{64}$/.test(value);
export function validateSavedIndex(index,data,hash){
  if(index?.data_sha256!==hash)throw new Error('정밀 열분포와 형상 파일 버전 불일치');
  const expected=Object.entries(data.baseline.profiles).map(([source,profile])=>({source,profile,identity:null}));
  for(const series of data.series)if(series.final)expected.push({source:series.source,profile:series.final.profile,identity:series.final.id});
  for(const row of expected){
    const entry=index.profiles?.[row.profile];
    if(!entry||entry.source!==row.source||entry.identity!==row.identity||!validSHA(entry.sha256)||entry.profile!==`full-saved-profiles/${row.profile.split('/').pop()}`)throw new Error('기준/정밀 열분포 누락 또는 후보 불일치');
  }
  return index;
}
