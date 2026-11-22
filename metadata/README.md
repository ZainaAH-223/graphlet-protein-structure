# metadata

Small text files that describe the data. Committed to Git.

| File | What it is |
|---|---|
| `sources.tsv` | Every raw file: where it comes from, version, license. Read by `scripts/download.py`. |
| `checksums.sha256` | sha256 of each raw file. Written by `scripts/download.py` on the first download, checked on every later run. Verify by hand with `shasum -a 256 -c metadata/checksums.sha256`. |

Planned (added when created): `domains.tsv` (CATH domain set), `idp_proteins.tsv` (disorder extension).

## Data dictionary

### sources.tsv
Tab-separated, one row per raw file.

| Column | Meaning | Example |
|---|---|---|
| `file` | Path where the file is saved, relative to the project folder | `data/raw/1MBN.pdb` |
| `source` | Database or provider | `RCSB PDB` |
| `url` | Exact download URL | `https://files.rcsb.org/download/1MBN.pdb` |
| `version` | Database release, or `-` if the entry is not versioned | `-` |
| `license` | License of the data | `CC0` |

### checksums.sha256
Two columns separated by two spaces (standard `shasum` format): the sha256 hash, then the file path.
