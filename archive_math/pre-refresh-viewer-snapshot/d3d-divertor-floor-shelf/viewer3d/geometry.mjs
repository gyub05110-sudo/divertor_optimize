export function revolve(r,z,phi){return [r*Math.cos(phi),z,r*Math.sin(phi)];}
export function heatAt(samples,s){
 if(!samples.length||s<samples[0][0]||s>samples.at(-1)[0])return null;
 let lo=0,hi=samples.length-1;
 while(lo<hi){const mid=(lo+hi)>>1;if(samples[mid][0]<s)lo=mid+1;else hi=mid;}
 if(Math.abs(samples[lo][0]-s)<1e-10)return samples[lo][1];
 if(lo===0)return null;
 const a=samples[lo-1],b=samples[lo];
 if(b[0]-a[0]>.02)return null;
 return a[1]+(b[1]-a[1])*(s-a[0])/(b[0]-a[0]);
}
export function densify(points,step=.005){
 const out=[[...points[0],0]];let s=0;
 for(let i=1;i<points.length;i++){
  const a=points[i-1],b=points[i],d=Math.hypot(b[0]-a[0],b[1]-a[1]),n=Math.max(1,Math.ceil(d/step));
  for(let j=1;j<=n;j++)out.push([a[0]+(b[0]-a[0])*j/n,a[1]+(b[1]-a[1])*j/n,s+d*j/n]);
  s+=d;
 }return out;
}
