export function revolve(r,z,phi){return [r*Math.cos(phi),z,r*Math.sin(phi)];}
export function heatAt(samples,s){
 if(!samples.length||s<samples[0][0]||s>samples.at(-1)[0])return null;
 let lo=0,hi=samples.length-1;
 while(lo<hi){const mid=(lo+hi)>>1;if(samples[mid][0]<s)lo=mid+1;else hi=mid;}
 if(Math.abs(samples[lo][0]-s)<1e-10)return samples[lo][1];
 if(lo===0)return null;
 const a=samples[lo-1],b=samples[lo];
 if(b[0]-a[0]>.02||a[2]!==b[2])return null;
 return a[1]+(b[1]-a[1])*(s-a[0])/(b[0]-a[0]);
}
export function heatSpan(samples,a,b){
 const qa=heatAt(samples,a),qb=heatAt(samples,b);
 return qa===null||qb===null||heatAt(samples,(a+b)/2)===null?null:[qa,qb];
}
export function heatAtTriangle(vertices,values,point){
 const [a,b,c]=vertices,sub=(v,w)=>v.map((x,i)=>x-w[i]),dot=(v,w)=>v.reduce((sum,x,i)=>sum+x*w[i],0);
 const ab=sub(b,a),ac=sub(c,a),ap=sub(point,a),d00=dot(ab,ab),d01=dot(ab,ac),d11=dot(ac,ac),d20=dot(ap,ab),d21=dot(ap,ac);
 const denominator=d00*d11-d01*d01;if(denominator<=0)return null;
 const v=(d11*d20-d01*d21)/denominator,w=(d00*d21-d01*d20)/denominator,weights=[1-v-w,v,w];
 const sum=items=>items.reduce((total,item,i)=>total+weights[i]*item,0);
 return {q:values.every(Number.isFinite)?sum(values):null,R:sum(vertices.map(p=>Math.hypot(p[0],p[2]))),Z:sum(vertices.map(p=>p[1]))};
}
export function densify(points,step=.005){
 const out=[[...points[0],0]];let s=0;
 for(let i=1;i<points.length;i++){
  const a=points[i-1],b=points[i],d=Math.hypot(b[0]-a[0],b[1]-a[1]),n=Math.max(1,Math.ceil(d/step));
  for(let j=1;j<=n;j++)out.push([a[0]+(b[0]-a[0])*j/n,a[1]+(b[1]-a[1])*j/n,s+d*j/n]);
  s+=d;
}return out;
}

export function densifyWithHeat(points,samples,step=.005){
 const lengths=[0];for(let i=1;i<points.length;i++)lengths.push(lengths.at(-1)+Math.hypot(points[i][0]-points[i-1][0],points[i][1]-points[i-1][1]));
 const total=lengths.at(-1),base=densify(points,step);
 const anchors=[...new Set([...base.map(p=>p[2]),...samples.map(p=>p[0]).filter(s=>s>=0&&s<=total)])].sort((a,b)=>a-b);
 let edge=1;
 return anchors.map(s=>{
  while(edge<lengths.length-1&&s>lengths[edge])edge++;
  const a=points[edge-1],b=points[edge],span=lengths[edge]-lengths[edge-1],t=span?(s-lengths[edge-1])/span:0;
  return [a[0]+t*(b[0]-a[0]),a[1]+t*(b[1]-a[1]),s];
 });
}
