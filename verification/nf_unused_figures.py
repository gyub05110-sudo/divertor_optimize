"""Plot the observed, time-admitted endpoints of the held-out campaign."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from d3d_optimization.broadening_study import load_json


def plot(study):
    report=load_json(study/"summary.json")
    cases=list(dict.fromkeys(row["case_id"] for row in report["rows"]))
    fig,axes=plt.subplots(2,len(cases),figsize=(12.8,6.5),squeeze=False)
    colours={"direct":"#7c909c","screen":"#ca642b"}
    labels={"direct":"Direct fine search","screen":"Coarse selection + fine check"}
    for column,case in enumerate(cases):
        rows=[row for row in report["rows"] if row["case_id"]==case]
        x=np.arange(len(rows))
        axes[0,column].set_title(case,fontsize=12)
        for i,metric in enumerate(("peak_reduction_pct","A2_change_pct")):
            ax=axes[i,column]
            for j,policy in enumerate(("direct","screen")):
                values=[row[policy][metric] for row in rows]
                bars=ax.bar(x+(j-.5)*.36,values,width=.34,color=colours[policy],label=labels[policy])
                ax.bar_label(bars,labels=[f"{v:.2f}" for v in values],padding=3,fontsize=8)
            ax.axhline(0,color="#384750",lw=.8)
            ax.set_xticks(x,[str(row["seed"]) for row in rows],fontsize=9)
            ax.set_xlabel("Random seed",fontsize=10)
            ax.grid(axis="y",alpha=.18)
            ax.set_axisbelow(True)
            ax.spines[["top","right"]].set_visible(False)
            ax.margins(y=.25)
        axes[0,column].set_ylim(0,62)
    axes[0,0].set_ylabel("Peak reduction from baseline (%)")
    axes[1,0].set_ylabel("A₂ change from baseline (%)")
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",ncol=2,frameon=False,bbox_to_anchor=(.5,.94))
    fig.suptitle("200 s per method, including preparation",fontsize=15,y=.995)
    fig.text(.5,.015,"OMP 1 MW · fixed DIII-D wall · cubic1280 endpoints completed before the deadline",
             ha="center",fontsize=10)
    fig.subplots_adjust(top=.83,bottom=.12,left=.075,right=.99,wspace=.28,hspace=.5)
    for suffix in ("png","pdf"):
        fig.savefig(study/f"paired-results.{suffix}",dpi=180)
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--study",required=True,type=Path)
    args=p.parse_args()
    plot(args.study)

