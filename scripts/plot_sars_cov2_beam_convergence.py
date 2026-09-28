#!/usr/bin/env python3
"""Plot SARS-CoV-2 genome ensemble free energy vs. beam width.

Methodological convergence check, analogous to LinearCapR's Figure 7
(Otagaki2026-ev): does the beam-pruned ensemble free energy stabilize as
beam width increases, on a genome-scale sequence? Makes no claim about
which genomic positions are accessible or biologically important.

Reads benchmark/results/sars_cov2_beam_convergence_capb/out_beam*.txt (caps tied to the beam)
(each file's "# ensemble_free_energy_kcal_mol=..." comment line) and writes
manuscript/figures/figure_sars_cov2_beam_convergence.pdf.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = REPO_ROOT / "benchmark/results/sars_cov2_beam_convergence_capb"
OUTPUT = REPO_ROOT / "manuscript/figures/figure_sars_cov2_beam_convergence.pdf"

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"]
plt.rcParams["font.size"] = 14


def main() -> int:
    beams = []
    energies = []
    for path in sorted(RESULTS_DIR.glob("out_beam*.txt")):
        beam = int(re.search(r"out_beam(\d+)\.txt", path.name).group(1))
        text = path.read_text()
        m = re.search(r"# ensemble_free_energy_kcal_mol=(-?[\d.]+)", text)
        if not m:
            continue
        beams.append(beam)
        energies.append(float(m.group(1)))

    order = sorted(range(len(beams)), key=lambda i: beams[i])
    beams = [beams[i] for i in order]
    energies = [energies[i] for i in order]

    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.plot(beams, energies, "o-", color="#d62728", linewidth=2, markersize=8)
    ax.set_xscale("log")
    ax.set_xlabel("Beam width $b$", fontsize=16)
    ax.set_ylabel(r"Ensemble free energy $G_{\mathrm{ensemble}}$ (kcal/mol)", fontsize=16)
    ax.set_title(
        "Beam-width sensitivity of the SARS-CoV-2 genome ensemble\n(NC_045512.2, 29,903 nt, Turner 1999)",
        fontsize=14,
    )
    ax.set_xticks(beams)
    ax.set_xticklabels([str(b) for b in beams])
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUTPUT)
    print(f"Wrote {OUTPUT}")
    for b, e in zip(beams, energies):
        print(f"  beam={b}: {e:.2f} kcal/mol")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
