export function casesForFile(cases,sha256){return cases.filter(c=>c.shot_sha256===sha256);}
export async function fileSHA256(file){
 const digest=await crypto.subtle.digest('SHA-256',await file.arrayBuffer());
 return [...new Uint8Array(digest)].map(b=>b.toString(16).padStart(2,'0')).join('');
}
