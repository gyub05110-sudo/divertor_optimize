import sys,json,hashlib
from pathlib import Path
import numpy as np
from scipy.interpolate import BSpline
from shapely.geometry import Point,LineString,Polygon
from shapely.ops import nearest_points
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(r'C:/Users/FPL-user/Desktop/LAB/NEW AutoFRD/freegs_divertor_no-transport')
sys.path.insert(0,str(ROOT))
OUT=Path(__file__).parent
STUDY=Path('D:/DivertorOptimization/studies/regional-20261028/qualified/OMP')
m=json.loads((STUDY/'manifest.json').read_text())
surface=m['metadata']['surfaces'][3];a,b=np.array(surface['endpoints_RZ_m']);v=b-a;L=np.linalg.norm(v)
n=np.array([-v[1],v[0]])/L;tip=np.array([1.37256,-1.3247]);endR=1.4277
coeff=np.array(m['candidates']['selected']);knots=np.r_[np.zeros(4),np.linspace(0,1,6)[1:-1],np.ones(4)]

def floor(c,count=10001):
    u=np.linspace(0,1,count);sp=BSpline(knots,np.r_[0,0,c[12:16],0,0],3)
    return a+u[:,None]*v+sp(u)[:,None]*n,sp

def metrics(c):
    f,sp=floor(c);near=np.array(nearest_points(Point(tip),LineString(f))[1].coords[0])
    radii={}
    for u in [0.,1.]:
        curvature=L*abs(sp(u,2))/(L*L+sp(u,1)**2)**1.5
        radii[str(u)]=float(1/curvature) if curvature else None
    rr=np.linspace(tip[0],endR,2001);zz=np.interp(rr,f[::-1,0],f[::-1,1]);gap=tip[1]-zz
    left=np.array([1.0882576947643479,-1.2853220897608946])-b
    outward=a-b
    exact_opening=np.degrees(np.arccos(np.dot(left,outward)/np.linalg.norm(left)/np.linalg.norm(outward)))
    arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(f[::-1],axis=0),axis=1))]
    apparent={}
    for distance in [.005,.01,.02]:
        point=np.array([np.interp(distance,arc,f[::-1,j]) for j in [0,1]])
        ray=point-b
        apparent[str(distance)]=float(np.degrees(np.arccos(np.dot(left,ray)/np.linalg.norm(left)/np.linalg.norm(ray))))
    return dict(tip_floor_shortest_gap_mm=float(np.linalg.norm(tip-near)*1000),nearest_floor_RZ_m=near.tolist(),
        mouth_vertical_gap_mm=float(gap[0]*1000),cavity_cross_section_m2=float(np.trapezoid(gap,rr)),
        cavity_axisymmetric_volume_m3=float(np.trapezoid(2*np.pi*rr*gap,rr)),
        endpoint_displacement_m=sp([0,1]).tolist(),endpoint_displacement_slope=sp([0,1],1).tolist(),
        left_exact_opening_angle_deg=float(exact_opening),left_tangent_turn_deg=float(180-exact_opening),
        left_chord_opening_deg=apparent,floor_one_sided_endpoint_radius_m=radii),rr,gap

record={'scope':'Post-analysis only; optimized study and production model unchanged','reference':m['reference'],
    'input_sha256':{x:hashlib.sha256((STUDY/x).read_bytes()).hexdigest() for x in ['manifest.json','baseline_refined.npz','selected_refined.npz']},
    'geometry':{},'counterfactuals':{}}
fig,axes=plt.subplots(1,3,figsize=(16,5),layout='constrained')
for name,c,color,style in [('baseline',np.zeros(24),'#2876b5','--'),('selected',coeff,'#d66a27','-')]:
    d=np.load(STUDY/(name+'_refined.npz'));w=d['wall'];mtr,rr,gap=metrics(c)
    poly=Polygon(w);near=np.array(mtr['nearest_floor_RZ_m'])
    mtr['gap_midpoint_inside_wall']=bool(poly.contains(Point((tip+near)/2)))
    coarse=LineString(w[surface['segment_index']+3*639:surface['segment_index']+3*639+641])
    mtr['gap_641_vs_10001_nodes_difference_mm']=abs(coarse.distance(Point(tip))*1000-mtr['tip_floor_shortest_gap_mm'])
    q=d['q_MW_m2'];points=d['points'];mtr['peak_RZ_m']=points[np.argmax(q)].tolist()
    for label,center in [('left_corner',b),('fixed_tip',tip),('floor_below_tip',near)]:
        mask=np.linalg.norm(points-center,axis=1)<.005
        mtr[label+'_qmax_within_5mm_MW_m2']=float(q[mask].max()) if mask.any() else None
    record['geometry'][name]=mtr
    for ax in axes[:2]:ax.plot(w[:,0],w[:,1],style,color=color,lw=2,label=name)
    axes[1].plot([tip[0],near[0]],[tip[1],near[1]],color=color,lw=1.5,marker='o',ms=4)
    axes[1].annotate(f"{mtr['tip_floor_shortest_gap_mm']:.2f} mm",(near+tip)/2,xytext=(8,0),textcoords='offset points',color=color,fontsize=10)
    axes[2].plot(rr,gap*1000,style,color=color,label=name)
    np.savetxt(OUT/(name+'_gap_profile.csv'),np.c_[rr,gap*1000],delimiter=',',header='R_m,vertical_gap_mm',comments='')
axes[0].plot(*b,'ko',ms=4);axes[0].set(xlim=(1.15,1.205),ylim=(-1.365,-1.30),xlabel='R [m]',ylabel='Z [m]',title='Left junction: inherited corner + tight bend')
axes[1].set(xlim=(1.35,1.44),ylim=(-1.365,-1.305),xlabel='R [m]',ylabel='Z [m]',title='Right recess: fixed lip to moving floor')
for ax in axes[:2]:ax.set_aspect('equal')
axes[2].set(xlabel='R beneath fixed lip [m]',ylabel='Vertical gap [mm]',title='Recess opening profile',ylim=(0,40))
for ax in axes:ax.legend();ax.grid(alpha=.2)
fig.savefig(OUT/'corner_gap_analysis.png',dpi=180);plt.close(fig)
(OUT/'analysis.json').write_text(json.dumps(record,indent=2))
print('GEOMETRY',json.dumps(record['geometry'],indent=2),flush=True)

from d3d_optimization.boundary_evaluator import prepare_boundary
from d3d_optimization.reference_transport import ReferencePoint
e=prepare_boundary(ROOT,'OMP',ReferencePoint(**m['reference']),count=320,nodes=161,step=.0005,wall_step=.001,quadrature=64,bound=.1)
cases={'baseline':np.zeros(24),'selected':coeff.copy()}
cases['right_end_coefficient_zero']=coeff.copy();cases['right_end_coefficient_zero'][12]=0.
cases['left_end_coefficient_half']=coeff.copy();cases['left_end_coefficient_half'][15]*=.5
cases['both_local_changes']=cases['right_end_coefficient_zero'].copy();cases['both_local_changes'][15]*=.5
for label,c in cases.items():
    result=e.evaluate(c);row=e.describe(result);geom,_,_=metrics(c)
    record['counterfactuals'][label]={'coefficients_m':c.tolist(),'result':row,'geometry':geom,
        'scope':'320 donor CPU evaluation; no reoptimization or full four-level verification'}
    print('CASE',label,'feasible',row['feasible'],'peak',row.get('peak_MW_m2'),'gap',geom['tip_floor_shortest_gap_mm'],flush=True)
    (OUT/'analysis.json').write_text(json.dumps(record,indent=2))
print('DONE',flush=True)
