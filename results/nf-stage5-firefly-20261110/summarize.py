"""Publish whole-wall FIREFLY paired intervals and their exact tabular sources."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np


CASES = ("3264fdd71f", "7ce879f1b3", "70481becc0", "5d7f95a730", "1255ab8051", "d880609264", "f317e77cf7")
FORCED = ("floor-only", "floor-entry", "floor-shelf", "floor-connector")
BLUE, ORANGE, INK, GRAY = "#215d7c", "#a54a18", "#273039", "#67737a"
KOREAN = {"floor-only":"바닥", "floor-shelf":"바닥+선반", "floor-entry":"바닥+입구", "floor-connector":"바닥+연결부",
          "floor":"바닥", "shelf":"선반", "connector":"연결부", "entry":"입구", "fixed":"고정벽",
          "supported":"지지", "contradicted":"불일치", "unresolved":"유보"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def columns(pair, prefix):
    ff = pair["firefly"]
    lo,hi = ff["ci95_pct"]
    val = ff["estimate_pct"]
    require(all(math.isfinite(x) for x in (lo,hi,val)) and lo<=hi, "invalid paired interval")
    values = dict(pct=val,ci95_low_pct=lo,ci95_high_pct=hi,uncertainty_pp=ff["uncertainty_pp"],
                  direction=ff["direction"],ci_only_1pct=ff["ci_only_1pct"],precision_1pct=ff["precision_1pct"],
                  verdict_1pct=ff["verdict_1pct"],model_pct=pair["model_pct"],stage4_pct=pair["stage4_pct"])
    return {prefix+"_"+k:v for k,v in values.items()}


def tables(metrics, digest):
    require(metrics["schema"]=="nf-stage5-whole-wall-metrics.v1" and metrics["eligible_files"]==7 and
            metrics["original_registration_files"]==24 and metrics["role_count"]==51 and metrics["distinct_case_walls"]==46
            and metrics["ordered_bootstrap_draws"]==6**6 and metrics["unique_bootstrap_counts"]==462,
            "requires finalized, complete Stage 5 metrics")
    require(tuple(c["case_id"] for c in metrics["cases"])==CASES, "wrong case inventory")
    forecast,products = [],[]
    for case in metrics["cases"]:
        cid = case["case_id"]
        roles = {(r["seed"],r["role"]):r for r in case["roles"]}
        baseline = roles[30713,"baseline"]
        family = case["predicted_family"]
        require(family in FORCED, "unregistered family")
        seeds = [30713,30757] if cid=="70481becc0" else [30713]
        require([c["optimizer_seed"] for c in case["comparisons"]]==seeds, "missing or repeated product seed")
        require(len(roles)==len(case["roles"])==5+2*len(seeds), "incomplete role inventory")
        p = case["comparisons"][0]
        selected = p["predicted_forced_vs_best_forced"]
        require(selected["target"]==family and set(selected["references"])==set(FORCED),"incorrect forced comparator")
        r = dict(case_id=cid,optimizer_seed=30713,predicted_family=family,
                 best_forced_families="|".join(selected["firefly"]["best_references"]),
                 predicted_forced_wall_id=roles[30713,family]["wall_id"],floor_only_wall_id=roles[30713,"floor-only"]["wall_id"])
        for model in ("stage4","model","firefly"):
            base = baseline[model]
            require(base["peak_surface"]==case["baseline_hotspot"][model],"reported baseline hotspot differs")
            r[model+"_baseline_region"] = base["peak_surface"]
            r[model+"_baseline_peak_MW_m2"] = base["peak_MW_m2"]
        r["firefly_baseline_peak_ci95_low_MW_m2"],r["firefly_baseline_peak_ci95_high_MW_m2"] = baseline["firefly"]["peak_ci95_MW_m2"]
        for arm in FORCED:
            for model in ("stage4","model","firefly"):
                r[model+"_"+arm.replace("-","_")+"_peak_MW_m2"] = roles[30713,arm][model]["peak_MW_m2"]
        r.update(columns(selected,"predicted_vs_best") | columns(p["floor_only_vs_best_forced"],"floor_vs_best") |
                 columns(p["predicted_forced_vs_floor_only"],"predicted_vs_floor"),metrics_sha256=digest)
        forecast.append(r)
        for c in case["comparisons"]:
            seed = c["optimizer_seed"]
            row = dict(case_id=cid,optimizer_seed=seed)
            for arm in ("baseline","hotspot","full"):
                role = baseline if arm=="baseline" else roles[seed,arm]
                row[arm+"_wall_id"] = role["wall_id"]
                for model in ("stage4","model","firefly"):
                    row[model+"_"+arm+"_peak_MW_m2"] = role[model]["peak_MW_m2"]
            row.update(columns(c["hotspot_vs_full"],"auto_vs_full") | columns(c["hotspot_vs_baseline"],"auto_vs_baseline") |
                       columns(c["full_vs_baseline"],"full_vs_baseline"),metrics_sha256=digest)
            products.append(row)
    return forecast,products


def markdown(forecast, products):
    def val(row,metric):
        q,lo,hi = [row[metric+"_"+k] for k in ("pct","ci95_low_pct","ci95_high_pct")]
        return f"{q:+.2f} [{lo:+.2f}, {hi:+.2f}]"
    out = ["상대 피크 차이 (%)와 짝지은 Monte Carlo 95% 구간. 음수는 분자 형상의 피크가 더 낮음을 뜻한다.","",
           "| 평형 | 지목한 가족 | 지목/4종 최저 | 1% 판정 | 바닥/4종 최저 | 1% 판정 | 지목/바닥 |",
           "|---|---|---:|---|---:|---|---:|"]
    for r in forecast:
        direct = val(r,"predicted_vs_floor") if r["predicted_family"]!="floor-only" else "동일 형상"
        out.append(f'| {r["case_id"]} | {KOREAN[r["predicted_family"]]} | {val(r,"predicted_vs_best")} | '
                   f'{KOREAN[r["predicted_vs_best_verdict_1pct"]]} | {val(r,"floor_vs_best")} | '
                   f'{KOREAN[r["floor_vs_best_verdict_1pct"]]} | {direct} |')
    out.extend(["","| 평형 | 최적화 시드 | 자동/전체 (%) | 자동/원벽 (%) | 전체/원벽 (%) |",
                "|---|---:|---:|---:|---:|"])
    for r in products:
        out.append(f'| {r["case_id"]} | {r["optimizer_seed"]} | {val(r,"auto_vs_full")} | '
                   f'{val(r,"auto_vs_baseline")} | {val(r,"full_vs_baseline")} |')
    out.extend(["","| 원벽 평형 | 원래 모델 정밀 열점 | 2 mm 모델 열점·피크 (MW/m²) | FIREFLY 열점·피크 [95% 구간] (MW/m²) |",
                "|---|---|---:|---:|"])
    for r in forecast:
        surfaces = [KOREAN.get(r[m+"_baseline_region"],r[m+"_baseline_region"]) for m in ("stage4","model","firefly")]
        out.append(f'| {r["case_id"]} | {surfaces[0]} | {surfaces[1]} {r["model_baseline_peak_MW_m2"]:.3f} | '
                   f'{surfaces[2]} {r["firefly_baseline_peak_MW_m2"]:.3f} '
                   f'[{r["firefly_baseline_peak_ci95_low_MW_m2"]:.3f}, {r["firefly_baseline_peak_ci95_high_MW_m2"]:.3f}] |')
    return "\n".join(out)


def intervals(ax, rows, field, offset, color, marker, open_mark=False):
    y = np.arange(len(rows))+offset
    estimates = [r[field+"_pct"] for r in rows]
    lows,highs = ([r[field+key] for r in rows] for key in ("_ci95_low_pct","_ci95_high_pct"))
    ax.hlines(y,lows,highs,color=color,linewidth=1.3)
    ax.plot(estimates,y,linestyle="none",marker=marker,ms=5.3,color=color,
            markerfacecolor="white" if open_mark else color,markeredgewidth=1.15)


def chart(forecast, products):
    fig,axes = plt.subplots(2,2,figsize=(12.,8.7))
    fig.subplots_adjust(left=.132,right=.96,top=.842,bottom=.108,hspace=.62,wspace=.47)
    fig.suptitle("Whole-wall heat peaks on the fixed divertor designs",ha="left",x=.042,y=.978,fontsize=17,weight="semibold")
    fig.text(.042,.925,"7 eligible equilibria at 1400 ms  ·  1 MW outboard-midplane input  ·  matched 2 mm wall rings",
             fontsize=10,color=GRAY)
    first,base,full,wall = axes.ravel()
    labels = [r["case_id"] for r in forecast]
    product_labels = [r["case_id"]+("  *" if r["optimizer_seed"]==30757 else "") for r in products]
    for ax,title,ids in zip(axes.ravel(),("A   Forced family vs lowest of four","B   Original wall",
                        "C   Automatic vs full search","D   Automatic vs original wall"),
                        (labels,labels,product_labels,product_labels)):
        ax.set_title(title,loc="left",fontsize=11,weight="semibold",pad=24)
        ax.set_yticks(np.arange(len(ids)),ids,fontfamily="DejaVu Sans Mono",fontsize=8.4)
        ax.set_ylim(len(ids)-.43,-.56)
        ax.spines[["top","right","left"]].set_visible(False)
        ax.spines["bottom"].set_color("#a5afb5")
        ax.tick_params(axis="y",length=0,pad=6)
        ax.tick_params(axis="x",color="#a5afb5",labelsize=8.5)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.grid(axis="x",linewidth=.55,color="#dce2e5")
        ax.set_axisbelow(True)
        ax.margins(x=.07)
    intervals(first,forecast,"predicted_vs_best",-.16,BLUE,"o")
    intervals(first,forecast,"floor_vs_best",.16,ORANGE,"s",True)
    first.axvline(1.,color=INK,ls=(0,(4,3)),lw=.95,zorder=1)
    first.set_xlabel("Peak above lowest forced design (%)")
    def mark(color,marker,label,open_mark=False):
        return Line2D([0],[0],ls="none",marker=marker,color=color,markerfacecolor="white" if open_mark else color,
                      label=label,ms=5.2)
    def legend(ax,handles):
        ax.legend(handles=handles,loc="lower left",bbox_to_anchor=(-.008,1.013),ncol=3,frameon=False,fontsize=8,
                  handlelength=1.15,handletextpad=.35,columnspacing=.95,borderpad=0,borderaxespad=0)
    legend(first,[mark(BLUE,"o","Predicted family"),mark(ORANGE,"s","Always floor",True),
                  Line2D([0],[0],ls="--",lw=1.,color=INK,label="1% limit")])
    ys = np.arange(len(forecast))
    base.hlines(ys+.15,[r["firefly_baseline_peak_ci95_low_MW_m2"] for r in forecast],
                [r["firefly_baseline_peak_ci95_high_MW_m2"] for r in forecast],lw=1.3,color=BLUE)
    base.plot([r["firefly_baseline_peak_MW_m2"] for r in forecast],ys+.15,"o",color=BLUE,ms=5.3)
    base.plot([r["model_baseline_peak_MW_m2"] for r in forecast],ys-.15,"D",color=ORANGE,mfc="white",ms=5.1)
    base.set_xlabel("Whole-wall peak (MW/m²)")
    legend(base,[mark(BLUE,"o","FIREFLY"),mark(ORANGE,"D","Reduced model",True)])
    for ax,name in ((full,"auto_vs_full"),(wall,"auto_vs_baseline")):
        intervals(ax,products,name,0.,BLUE,"o")
        ax.axvline(0.,ls=(0,(4,3)),color=INK,lw=.9,zorder=1)
        ax.set_xlabel("Automatic peak relative to reference (%)")
        ax.text(0,1.025,"Negative values favor automatic search",transform=ax.transAxes,fontsize=8,color=GRAY)
    fig.text(.042,.045,"Points show pooled peak estimates. Bars show 95% paired seed-bootstrap intervals. "
             "* Extra automatic/full pair from optimizer seed 30757. Other designs use seed 30713.",fontsize=8.7,color=GRAY)
    return fig


def publish(study, output):
    data = (study/"metrics.json").read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    forecast,products = tables(json.loads(data),digest)
    output.mkdir(parents=True,exist_ok=True)
    for name,rows in (("firefly-forecast",forecast),("firefly-products",products)):
        with (output/(name+".csv")).open("w",encoding="utf-8",newline="") as file:
            writer=csv.DictWriter(file,fieldnames=rows[0].keys(),lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    with plt.rc_context({"font.family":"DejaVu Sans","font.size":9.3,"pdf.fonttype":42,"axes.labelcolor":INK,
                         "text.color":INK,"xtick.color":INK,"ytick.color":INK,"figure.facecolor":"white"}):
        fig = chart(forecast,products)
        try:
            subject = "Stage 5 whole-wall comparison. Source metrics.json SHA-256 "+digest
            fig.savefig(output/"firefly-comparisons.pdf",metadata={"Subject":subject,"CreationDate":None,"ModDate":None})
            fig.savefig(output/"firefly-comparisons.png",dpi=300,metadata={"Description":subject})
        finally:
            plt.close(fig)
    return forecast,products


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--study",type=Path,required=True)
    p.add_argument("--output",type=Path,default=Path(__file__).resolve().parent)
    a=p.parse_args()
    print(markdown(*publish(a.study,a.output)))
