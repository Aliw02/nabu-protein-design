import argparse,hashlib,json,re
from itertools import combinations
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr

def fnv(s):
 h=2166136261
 for b in s.encode(): h=((h^b)*16777619)&0xffffffff
 return h
def parse(s):
 z=re.findall(r"[A-Z][0-9]+[A-Z*]",str(s))
 return tuple(sorted(z,key=lambda x:(int(x[1:-1]),x)))
def pct(a,x):
 a=np.asarray(a,float); return float((np.sum(a<x)+.5*np.sum(a==x))/len(a))
def zscore(x):
 x=np.asarray(x,float); s=x.std()
 return np.zeros_like(x) if s==0 else (x-x.mean())/s
def main(path,out):
 o=Path(out);o.mkdir(parents=True,exist_ok=True)
 d=pd.read_csv(path)
 mc=next((c for c in ["mutant","mutation","mutations","variant"] if c in d.columns),None)
 yc=next((c for c in ["DMS_score","score","fitness","mean","Fitness"] if c in d.columns),None)
 if not mc or not yc: raise SystemExit("Need mutation and fitness columns; found "+str(list(d.columns)))
 d=d[[mc,yc]].dropna().copy(); d["ms"]=d[mc].map(parse)
 d=d[d.ms.map(lambda x:3<=len(x)<=5 and "*" not in "".join(x))].copy()
 d["cid"]=d.ms.map(lambda x:":".join(x));d["bucket"]=d.cid.map(lambda x:fnv(x)%10)
 v=d[d.bucket<=6].copy();h=d[d.bucket>=7].copy()
 g=float(v[yc].mean()); C={}
 for m in sorted({m for x in v.ms for m in x}):
  z=v[v.ms.map(lambda x:m in x)][yc];n=len(z);C[m]=(float(z.mean())-g)*n/(n+1)
 def additive(ms):return g+sum(C.get(m,0) for m in ms)
 v["b2"]=v.ms.map(additive);v["resid"]=v[yc]-v.b2;P={}
 for p in sorted({p for x in v.ms for p in combinations(x,2)}):
  z=v[v.ms.map(lambda x:set(p).issubset(x))].resid;n=len(z)
  if n:P[p]=(float(z.mean())*n/(n+1),n)
 def raw(ms):
  return additive(ms)+sum(P[p][0] for p in combinations(ms,2) if p in P)
 e=h[h.ms.map(lambda x:all(m in C for m in x) and any(p in P for p in combinations(x,2)))].copy()
 e["B2_ADDITIVE"]=e.ms.map(additive);e["B3_RAW_PAIR"]=e.ms.map(raw)
 # Frozen consensus: equal standardized additive + pair scores, with >=2-support pair evidence preferred.
 pair_supported=e.ms.map(lambda x:any(p in P and P[p][1]>=2 for p in combinations(x,2))).to_numpy()
 b5=(zscore(e.B2_ADDITIVE)+zscore(e.B3_RAW_PAIR))/2
 e["B5_NABU_CONSENSUS"]=np.where(pair_supported,b5,zscore(e.B2_ADDITIVE))
 # Deterministic visible-label rotation control; rebuild both memories.
 sv=v.copy();sv[yc]=np.roll(sv[yc].to_numpy(),161%len(sv));sg=float(sv[yc].mean());SC={}
 for m in sorted({m for x in sv.ms for m in x}):
  z=sv[sv.ms.map(lambda x:m in x)][yc];n=len(z);SC[m]=(float(z.mean())-sg)*n/(n+1)
 def sa(ms):return sg+sum(SC.get(m,0) for m in ms)
 sv["sb"]=sv.ms.map(sa);sv["sr"]=sv[yc]-sv.sb;SP={}
 for p in sorted({p for x in sv.ms for p in combinations(x,2)}):
  z=sv[sv.ms.map(lambda x:set(p).issubset(x))].sr;n=len(z)
  if n:SP[p]=(float(z.mean())*n/(n+1),n)
 e["SHUFFLED_B5"]=e.ms.map(lambda x:sa(x)+sum(SP[p][0] for p in combinations(x,2) if p in SP))
 # Phase 1: candidate freeze contains no hidden truth.
 cols=["B2_ADDITIVE","B3_RAW_PAIR","B5_NABU_CONSENSUS","SHUFFLED_B5"]
 freeze={"version":"NABU_V7_CANDIDATE_FREEZE","eligible_hidden":len(e),"arms":{}}
 for c in cols:
  q=e.sort_values([c,"cid"],ascending=[False,True]).head(50)
  freeze["arms"][c]=[{"candidate_id":r.cid,"predicted_score":float(getattr(r,c))} for r in q.itertuples()]
 rawj=json.dumps(freeze,sort_keys=True,separators=(",",":")).encode()
 sha=hashlib.sha256(rawj).hexdigest();freeze["sha256_without_self_hash"]=sha
 (o/"CANDIDATE_FREEZE_PRE_REVEAL.json").write_text(json.dumps(freeze,indent=2))
 (o/"CANDIDATE_FREEZE_SHA256.txt").write_text(sha+"\n")
 # Phase 2 reveal/evaluation.
 ys=e[yc].to_numpy(float);cut=np.quantile(ys,.99);best,worst=float(max(ys)),float(min(ys))
 def met(c):
  q=e.sort_values([c,"cid"],ascending=[False,True]);x=float(q.iloc[0][yc]);t5=q.head(5)[yc].to_numpy(float);t50=q.head(50)[yc].to_numpy(float)
  return {"spearman":float(spearmanr(e[c],e[yc]).statistic),"top1_percentile":pct(ys,x),"top5_mean_percentile":float(np.mean([pct(ys,a) for a in t5])),"top50_top1pct_enrichment":float(np.mean(t50>=cut)/.01),"normalized_top1_regret":float((best-x)/(best-worst)) if best>worst else 0}
 R={"version":"NABU_V7_RESULT","rows":len(d),"visible":len(v),"hidden":len(h),"eligible_hidden":len(e),"freeze_sha256":sha}
 for c in cols:R[c]=met(c)
 A,B=R["B5_NABU_CONSENSUS"],R["B3_RAW_PAIR"]
 wins=sum([A["top1_percentile"]>B["top1_percentile"],A["top5_mean_percentile"]>B["top5_mean_percentile"],A["top50_top1pct_enrichment"]>B["top50_top1pct_enrichment"],A["normalized_top1_regret"]<B["normalized_top1_regret"]])
 R["b5_vs_raw_pair_design_metric_wins"]=wins
 R["preregistered_pass"]=bool(wins>=3 and A["top50_top1pct_enrichment"]>1 and abs(R["SHUFFLED_B5"]["spearman"])<.05)
 (o/"V7_RESULTS.json").write_text(json.dumps(R,indent=2));e.drop(columns=["ms"]).to_csv(o/"V7_REVEALED_PREDICTIONS.csv",index=False)
 print(json.dumps(R,indent=2))
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("csv");p.add_argument("--out",default="nabu_v7_results");a=p.parse_args();main(a.csv,a.out)
