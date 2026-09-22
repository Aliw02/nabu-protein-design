import argparse,hashlib,json,re
from itertools import combinations
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
def H(s):
 h=2166136261
 for b in s.encode(): h=((h^b)*16777619)&0xffffffff
 return h
def M(s): return tuple(sorted(re.findall(r"[A-Z]\d+[A-Z]",str(s)),key=lambda x:int(x[1:-1])))
def pct(a,x): a=np.asarray(a); return float((sum(a<x)+.5*sum(a==x))/len(a))
a=argparse.ArgumentParser();a.add_argument("csv");a.add_argument("--out",default="results_v4");q=a.parse_args()
o=Path(q.out);o.mkdir(exist_ok=True)
d=pd.read_csv(q.csv); mc=next(c for c in ["mutant","mutation","mutations"] if c in d); yc=next(c for c in ["DMS_score","score","fitness","mean"] if c in d)
d=d[[mc,yc]].dropna();d["ms"]=d[mc].map(M);d=d[d["ms"].map(lambda x:2<=len(x)<=5)];d["b"]=d["ms"].map(lambda x:H(":".join(x))%10)
v=d[d.b<=6].copy(); h=d[d.b>=7].copy(); g=float(v[yc].mean()); C={}
for m in {m for x in v["ms"] for m in x}:
 z=v[v["ms"].map(lambda x:m in x)][yc];n=len(z);C[m]=(float(z.mean())-g)*n/(n+1)
base=lambda x:g+sum(C.get(m,0) for m in x)
v["base"]=v["ms"].map(base);v["r"]=v[yc]-v.base;P={}
for p in {p for x in v["ms"] for p in combinations(x,2)}:
 z=v[v["ms"].map(lambda x:set(p).issubset(x))].r;n=len(z)
 if n:P[p]=(float(z.mean())*n/(n+1),n/(n+1))
def full(x):
 z=[P[p] for p in combinations(x,2) if p in P]
 return base(x)+(sum(a*b for a,b in z)/sum(b for a,b in z) if z else 0)
e=h[h["ms"].map(lambda x:all(m in C for m in x))].copy();e["component_score"]=e["ms"].map(base);e["lightning_score"]=e["ms"].map(full)
f=e.sort_values("lightning_score",ascending=False)[[mc,"component_score","lightning_score"]].head(50);f.to_csv(o/"candidate_manifest_pre_reveal.csv",index=False)
sha=hashlib.sha256((o/"candidate_manifest_pre_reveal.csv").read_bytes()).hexdigest();(o/"candidate_manifest_sha256.txt").write_text(sha+"\n")
ys=e[yc].to_numpy();cut=np.quantile(ys,.95)
def met(c):
 z=e.sort_values(c,ascending=False);x=float(z.iloc[0][yc]);t5=z.head(5)[yc];t50=z.head(50)[yc];best=max(ys);worst=min(ys)
 return {"spearman":float(spearmanr(e[c],e[yc]).statistic),"top1_percentile":pct(ys,x),"top5_mean_percentile":float(np.mean([pct(ys,x) for x in t5])),"top50_top5pct_enrichment":float(np.mean(t50>=cut)/.05),"normalized_top1_regret":float((best-x)/(best-worst))}
R={"rows":len(d),"visible":len(v),"hidden":len(h),"eligible_hidden":len(e),"manifest_sha256":sha,"component_memory":met("component_score"),"full_lightning":met("lightning_score")}
sv=v.copy();sv[yc]=np.roll(sv[yc].to_numpy(),161%len(sv));sg=float(sv[yc].mean());SC={}
for m in {m for x in sv["ms"] for m in x}:
 z=sv[sv["ms"].map(lambda x:m in x)][yc];n=len(z);SC[m]=(float(z.mean())-sg)*n/(n+1)
e["shuffle"]=e["ms"].map(lambda x:sg+sum(SC.get(m,0) for m in x));R["shuffle_spearman"]=float(spearmanr(e["shuffle"],e[yc]).statistic)
L,B=R["full_lightning"],R["component_memory"];wins=sum([L["top1_percentile"]>B["top1_percentile"],L["top5_mean_percentile"]>B["top5_mean_percentile"],L["top50_top5pct_enrichment"]>B["top50_top5pct_enrichment"],L["normalized_top1_regret"]<B["normalized_top1_regret"]]);R["design_metric_wins"]=wins;R["preregistered_pass"]=bool(wins>=3 and L["top50_top5pct_enrichment"]>1 and abs(R["shuffle_spearman"])<.05)
(o/"results.json").write_text(json.dumps(R,indent=2));e.drop(columns=["ms"]).to_csv(o/"revealed_predictions.csv",index=False);print(json.dumps(R,indent=2))

