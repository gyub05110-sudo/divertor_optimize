export function revolve([r,z],phi){return [r*Math.cos(phi),z,r*Math.sin(phi)];}
export function segmentsForWall(wall,lowerOnly=true){
  const out=[];
  for(let i=0;i<wall.length-1;i++){
    const a=wall[i],b=wall[i+1];
    if(!lowerOnly || (a[1]<-1.05 && b[1]<-1.05))out.push([a,b]);
  }
  return out;
}
export function selectedRecord(series,index,final){return final?series.final??null:series.stages[index-1]??null;}
export function legLayers(data,record,display){
  if(!data)return [];
  const result=[],add=(id,ghost)=>{if(data.cases[id])result.push({legs:data.cases[id],ghost});};
  if(display==='baseline'||!record)add('baseline',false);
  else {if(display==='overlay')add('baseline',true);add(record.id,false);}
  return result;
}
export function maxTangentError(c){return Math.max(...(c.tangent_angle_errors_deg??[c.tangent_angle_error_deg]));}
export function minRadius(c){return Math.min(...Object.values(c.actual_minimum_radius_m).filter(v=>typeof v==='number'&&Number.isFinite(v)));}
export function heatRGB(q,max){
  const stops=[[.08,.11,.18],[.29,.13,.46],[.66,.20,.38],[.95,.42,.16],[1,.77,.22],[1,.97,.66]];
  const t=Math.max(0,Math.min(1,q/max))*(stops.length-1),i=Math.min(stops.length-2,Math.floor(t)),f=t-i;
  return stops[i].map((v,k)=>v*(1-f)+stops[i+1][k]*f);
}
