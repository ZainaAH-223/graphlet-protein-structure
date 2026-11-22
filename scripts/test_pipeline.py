"""Pipeline test: run the whole pipeline on one known protein.

1MBN -> Cα contact network -> plot -> ORCA orbit counts -> TM-align vs 2HHB chain A.
Run from the project folder:  pixi run test-pipeline
"""
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no screen needed (works on GitHub Actions too)
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import yaml
from Bio.PDB import PDBParser, is_aa
from Bio.PDB.Polypeptide import index_to_one, three_to_index
from tmtools import tm_align

cfg = yaml.safe_load(open("config/config.yaml"))
cutoff = cfg["network"]["cutoff"]
p0 = cfg["pipeline_test"]
OUT = Path("results/pipeline_test")
OUT.mkdir(parents=True, exist_ok=True)


def ca_coords(path, chain_id=None):
    """Cα coordinates + sequence of one chain (first chain by default), standard residues only."""
    model = next(iter(PDBParser(QUIET=True).get_structure("s", path)))
    chain = model[chain_id] if chain_id else next(iter(model))
    coords, seq = [], []
    for res in chain:
        if is_aa(res, standard=True) and "CA" in res:
            coords.append(res["CA"].coord)
            seq.append(index_to_one(three_to_index(res.get_resname())))
    return np.array(coords, dtype=float), "".join(seq)


def contact_network(coords, cutoff):
    """Nodes = residues; edge if Cα–Cα distance <= cutoff."""
    dist = np.linalg.norm(coords[:, None] - coords[None, :], axis=-1)
    i, j = np.where(np.triu(dist <= cutoff, k=1))
    G = nx.Graph()
    G.add_nodes_from(range(len(coords)))
    G.add_edges_from(zip(i.tolist(), j.tolist()))
    return G, dist


def run_orca(G, name):
    """ORCA node orbit counts -> array (n_nodes x 15)."""
    in_path, out_path = OUT / f"{name}.in", OUT / f"{name}.orca"
    with open(in_path, "w") as f:
        f.write(f"{G.number_of_nodes()} {G.number_of_edges()}\n")
        for u, v in sorted(G.edges()):
            f.write(f"{u} {v}\n")
    subprocess.run([cfg["graphlets"]["orca"], "node", str(cfg["graphlets"]["size"]),
                    str(in_path), str(out_path)], check=True, capture_output=True)
    return np.loadtxt(out_path, dtype=int, ndmin=2)


# 1. Contact network
coords, seq = ca_coords(p0["test_pdb"])
G, dist = contact_network(coords, cutoff)
degrees = np.array([d for _, d in sorted(G.degree())])
print(f"1MBN: {len(coords)} residues, {G.number_of_edges()} contacts at {cutoff} Å, "
      f"connected={nx.is_connected(G)}, mean degree={degrees.mean():.2f}")

# 2. Plot: contact map + network coloured by sequence position
fig, ax = plt.subplots(1, 2, figsize=(12, 5.5))
ax[0].imshow(dist <= cutoff, cmap="Greys", origin="lower")
ax[0].set(title=f"1MBN contact map (Cα, {cutoff} Å)", xlabel="residue", ylabel="residue")
pos = nx.spring_layout(G, seed=cfg["seed"])
nx.draw_networkx_edges(G, pos, ax=ax[1], alpha=0.2, width=0.5)
nodes = nx.draw_networkx_nodes(G, pos, ax=ax[1], node_size=15,
                               node_color=range(len(coords)), cmap="viridis")
fig.colorbar(nodes, ax=ax[1], label="residue index")
ax[1].set_title("1MBN residue contact network")
ax[1].axis("off")
fig.tight_layout()
fig.savefig(OUT / "1mbn_network.png", dpi=150)
print(f"plot -> {OUT / '1mbn_network.png'}")

# 3. ORCA: orbit 0 must equal each residue's degree
orbits = run_orca(G, "1mbn")
assert orbits.shape == (len(coords), 15), orbits.shape
assert (orbits[:, 0] == degrees).all(), "orbit 0 should equal degree"
print(f"ORCA: {orbits.shape[0]} residues x {orbits.shape[1]} orbits, orbit 0 == degree ✓")

# 4. TM-align: same fold should give TM-score > 0.5
coords2, seq2 = ca_coords(p0["compare_pdb"], p0["compare_chain"])
res = tm_align(coords, coords2, seq, seq2)
print(f"TM-align 1MBN vs 2HHB:A  TM-score={res.tm_norm_chain1:.3f}  RMSD={res.rmsd:.2f} Å")
assert res.tm_norm_chain1 > 0.5, "myoglobin and haemoglobin should share a fold"
print("Pipeline test passed.")
