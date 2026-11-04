"""Strict finite-list selection. Numerical qualification alone is not improvement."""
import math

EXPECTED_LEVELS=[('cubic',320),('cubic',640),('cubic',1280),('quintic',640)]


def ranked_improvements(rows):
    eligible=[]
    for row in rows:
        d=row['description'];peak=d.get('peak_MW_m2',float('nan'))
        if (d.get('feasible') is True and d.get('geometry_screen',{}).get('passed') is True
                and math.isfinite(peak) and peak>0 and 0<d.get('peak_ratio',float('inf'))<1):
            eligible.append(row)
    return sorted(eligible,key=lambda row:(row['description']['peak_MW_m2'],row['identity']))


def passes_improvement(v):
    levels=v.get('levels',[]);verdict=v.get('verdict',{});wall=v.get('wall_only_refinement',{})
    if (v.get('complete') is not True or verdict.get('resolution_screen_passed') is not True
            or verdict.get('improvement_confirmed') is not True
            or [(x.get('method'),x.get('tubes')) for x in levels]!=EXPECTED_LEVELS):return False
    for level in levels:
        geom=level.get('geometry_audit',{})
        if geom.get('passed_independent') is not True or geom.get('geometry_screen',{}).get('passed') is not True:return False
        for name in ('selected','baseline'):
            row=level.get('rows',{}).get(name,{})
            if row.get('audit_passed') is not True or row.get('feasible') is not True:return False
            peak=row.get('peak_MW_m2',float('nan'))
            if not math.isfinite(peak) or peak<=0:return False
    if wall.get('geometry_audit',{}).get('passed_independent') is not True:return False
    if not all(wall.get('rows',{}).get(k,{}).get('audit_passed') is True
               and wall['rows'][k].get('feasible') is True for k in ('selected','baseline')):return False
    fine=levels[2]['rows']
    return fine['selected']['peak_MW_m2']<fine['baseline']['peak_MW_m2']
