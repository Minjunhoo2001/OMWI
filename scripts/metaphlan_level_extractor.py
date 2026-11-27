#!/usr/bin/env python3
"""Extract the six taxonomic ranks from a merged MetaPhlAn abundance table."""
import argparse
from pathlib import Path
import pandas as pd

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("input", type=Path)
p.add_argument("output_dir", type=Path)
a = p.parse_args()
with a.input.open(encoding="utf-8-sig") as fh:
    header = None
    for i, line in enumerate(fh):
        first = line.split("\t", 1)[0].strip().lstrip("#").strip()
        if first in {"clade_name", "ID", "Taxon"}:
            header = i
            break
    if header is None:
        raise ValueError("No clade_name / ID / Taxon header found")
df = pd.read_csv(a.input, sep="\t", skiprows=header)
if df.shape[1] < 2:
    raise ValueError("Expected a taxon column and at least one sample column")
a.output_dir.mkdir(parents=True, exist_ok=True)
name = df.columns[0]
clades = df[name].astype(str)
for prefix, rank in zip("kpcofgs", ("kingdom", "phylum", "class", "order", "family", "genus", "species")):
    selected = df[clades.str.contains(r"(?:^|\|)" + prefix + r"__[^|]+$", regex=True)]
    if len(selected):
        selected.to_csv(a.output_dir / (a.input.stem + "_" + rank + ".tsv"), sep="\t", index=False)
