"""Dataset selection: filter CATH S40 domains, find folds (CATH topologies) that
meet GRAFENE's minimum category size, pick a balanced set of folds per class, and
write a ranked candidate list per fold.

Run from the project folder:  pixi run select-dataset
Writes: results/dataset/fold_counts.tsv   (all folds, with counts)
        metadata/domain_candidates.tsv    (ranked candidates per chosen fold; committed)
Why candidates and not the final 30: some domains will fail the structure checks
done later (X-ray only, chain breaks, network connectivity). The final set takes
the first `per_fold` candidates of each fold that pass, in rank order.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

full_cfg = yaml.safe_load(open("config/config.yaml"))
cfg = full_cfg["dataset"]
SEED = full_cfg["seed"]
MIN = cfg["min_per_category"]
OUT = Path("results/dataset")
OUT.mkdir(parents=True, exist_ok=True)

# CATH List File (CLF) columns, from README-cath-list-file-format.txt
CLF_COLUMNS = ["domain", "C", "A", "T", "H", "S35", "S60", "S95", "S100",
               "S100_count", "length", "resolution"]

domains = pd.read_csv(cfg["domain_list"], sep=r"\s+", comment="#",
                      header=None, names=CLF_COLUMNS)
s40 = set(Path(cfg["s40_list"]).read_text().split())
print(f"CATH domain list: {len(domains)} domains; S40 list: {len(s40)} domains")

# 1. Filters, applied in order, with counts so every step is visible
steps = [
    ("in S40 non-redundant set", domains["domain"].isin(s40)),
    (f"class in {cfg['classes']}", domains["C"].isin(cfg["classes"])),
    (f"length >= {cfg['min_length']}", domains["length"] >= cfg["min_length"]),
    (f"resolution <= {cfg['max_resolution']} Å (drops NMR = 999, obsolete = 1000)",
     domains["resolution"] <= cfg["max_resolution"]),
]
keep = pd.Series(True, index=domains.index)
print("\nFilters:")
for label, mask in steps:
    keep &= mask
    print(f"  after {label}: {keep.sum()}")
pool = domains[keep].copy()
pool["fold"] = pool[["C", "A", "T"]].astype(str).agg(".".join, axis=1)
pool["superfamily"] = pool[["C", "A", "T", "H"]].astype(str).agg(".".join, axis=1)

# 2. Domains per fold, with fold names (CATH names file: node, domain, ":name")
names = {}
for line in Path(cfg["names"]).read_text().splitlines():
    if line.startswith("#") or not line.strip():
        continue
    node, rest = line.split(None, 1)
    names[node] = rest.split(":", 1)[1].strip() if ":" in rest else ""
folds = (pool.groupby(["C", "fold"]).size().rename("domains").reset_index()
         .sort_values(["C", "domains", "fold"], ascending=[True, False, True]))
folds["name"] = folds["fold"].map(names).fillna("")
folds["qualifies"] = folds["domains"] >= MIN
folds.to_csv(OUT / "fold_counts.tsv", sep="\t", index=False)

# 3. Qualifying folds per class (GRAFENE: a category needs >= MIN networks)
q = folds[folds["qualifies"]]
per_class = q.groupby("C").size().reindex(cfg["classes"], fill_value=0)
print(f"\nFolds with >= {MIN} domains, per class:")
for c in cfg["classes"]:
    print(f"  class {c}: {per_class[c]} folds")
eligible = [c for c in cfg["classes"] if per_class[c] > 0]
if len(eligible) < 2:
    raise SystemExit("Fewer than 2 classes qualify: GRAFENE needs >= 2 categories per level.")
dropped = [c for c in cfg["classes"] if c not in eligible]
if dropped:
    print(f"  dropped by the >= {MIN} rule: class {dropped}")

# 4. Balanced design: k folds per class (auto = largest number every class can supply)
k_rule = int(per_class[eligible].min())
k = k_rule if cfg["folds_per_class"] == "auto" else int(cfg["folds_per_class"])
if k > k_rule:
    raise SystemExit(f"folds_per_class={k}, but only {k_rule} folds qualify in every class")
chosen = q[q["C"].isin(eligible)].groupby("C").head(k)   # most populated first
print(f"\nk = {k} folds per class ({'rule' if cfg['folds_per_class'] == 'auto' else 'set in config'}; "
      f"rule gives {k_rule}) -> {k * len(eligible)} folds x {cfg['per_fold']} = "
      f"{k * len(eligible) * cfg['per_fold']} domains")

# 5. Ranked candidates per fold: a seeded shuffle, one independent random stream per fold
rows = []
for _, f in chosen.iterrows():
    members = pool[pool["fold"] == f["fold"]].sort_values("domain")
    c, a, t = (int(x) for x in f["fold"].split("."))
    rng = np.random.default_rng([SEED, c, a, t])
    order = rng.permutation(len(members))
    ranked = members.iloc[order].head(cfg["candidates_per_fold"]).copy()
    ranked["rank"] = range(1, len(ranked) + 1)
    rows.append(ranked)
cand = pd.concat(rows)
cand["fold_name"] = cand["fold"].map(names).fillna("")
cand = cand[["domain", "C", "fold", "fold_name", "superfamily", "length",
             "resolution", "rank"]].rename(columns={"C": "class"})
Path("metadata").mkdir(exist_ok=True)
cand.to_csv("metadata/domain_candidates.tsv", sep="\t", index=False)

print(f"\nWrote metadata/domain_candidates.tsv: {len(cand)} candidates "
      f"(up to {cfg['candidates_per_fold']} per fold)")
print(cand.groupby("class").agg(folds=("fold", "nunique"), candidates=("domain", "size"),
                                median_length=("length", "median")).to_string())
