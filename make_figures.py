#!/usr/bin/env python3
"""Regenerate figures from archived experiment JSON, without refitting."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
HERE=Path(__file__).resolve().parent
p=json.loads((HERE/"paired_region_experiment_results.json").read_text())
d=json.loads((HERE/"homd_deep_results.json").read_text())
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"figure.dpi":150,"savefig.dpi":220})
fig,ax=plt.subplots(figsize=(5.5,3.6),layout="constrained")
values=[p["amp_hmt_accuracy"],p["full_hmt_accuracy"]]
bars=ax.bar(["V3–V4 fragment","Full-length 16S"],values,color=["#cc6677","#4477aa"],width=.6)
for bar in bars:ax.text(bar.get_x()+bar.get_width()/2,bar.get_height()+.015,f"{bar.get_height():.1%}",ha="center")
ci=p["bootstrap"]["accuracy_full_minus_amp_95pct_hmt_cluster_bootstrap"]
ax.text(.5,.37,f"Paired gain {p['absolute_accuracy_difference']:.1%}\n95% HMT-cluster interval {ci[0]:.1%}–{ci[2]:.1%}",ha="center",va="center",bbox={"facecolor":"white","edgecolor":"#999999"})
ax.set(ylabel="Correct HMT top-hit assignment",ylim=(0,1.02),title=f"Same {p['n_known']:,} known-HMT queries")
fig.savefig(HERE/"figure_1_paired_accuracy.png")
plt.close(fig)
random_auc=np.array([v["cos"]["known_detection_auroc"] for v in d["current_repeats"].values()])
temporal=d["temporal"]["cos"]["known_detection_auroc"]
fig,ax=plt.subplots(figsize=(5.2,3.5),layout="constrained")
ax.scatter(np.arange(1,6),random_auc,color="#4477aa",s=45,label="Within-v16 taxon holdout")
ax.hlines(temporal,.7,5.3,color="#cc6677",lw=2,label="v15 → v16")
ax.set(xticks=np.arange(1,6),xlabel="Within-release split",ylabel="Known versus new HMT AUROC",ylim=(.78,1.0),xlim=(.7,5.3))
ax.legend(loc="lower left",fontsize=8)
fig.savefig(HERE/"figure_2_split_sensitivity.png")
plt.close(fig)
