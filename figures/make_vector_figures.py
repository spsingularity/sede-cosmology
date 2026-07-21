"""Regenerate every paper figure as a vector PDF (alongside the existing PNGs).

matplotlib writes a true vector PDF when handed a ``*.pdf`` path, so this driver
monkeypatches ``Figure.savefig`` to emit a ``*.pdf`` next to every ``*.png`` it
writes, then runs the four figure-producing scripts as ``__main__``. The LaTeX
manuscript (paper/tex/main.tex) references the figures without an extension via
``\\graphicspath``, so pdflatex picks up the PDFs automatically.

Run from anywhere:
    python paper/make_vector_figures.py
"""
import os
import sys
import runpy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# src/ (for `import sede`) + src/experiments & src/scripts (flat-namespace sibling imports,
# e.g. run_e1_figure does `from run_lambda_verify import ...`).
for _p in (os.path.join(ROOT, "src", "scripts"), os.path.join(ROOT, "src", "experiments"), os.path.join(ROOT, "src")):
    sys.path.insert(0, _p)
os.chdir(ROOT)  # the run-scripts write to relative results/ paths

import matplotlib
matplotlib.use("Agg")
import matplotlib.figure as mfig

_orig_savefig = mfig.Figure.savefig


def _savefig(self, fname, *args, **kwargs):
    _orig_savefig(self, fname, *args, **kwargs)
    if isinstance(fname, (str, os.PathLike)):
        s = os.fspath(fname)
        if s.lower().endswith(".png"):
            pdf = s[:-4] + ".pdf"
            # dpi is irrelevant for vector output; keep bbox_inches etc.
            kw = {k: v for k, v in kwargs.items() if k != "dpi"}
            _orig_savefig(self, pdf, *args, **kw)
            print(f"  + vector  {pdf}")


mfig.Figure.savefig = _savefig

SCRIPTS = [
    "figures/figures.py",                         # fig1_mechanism .. fig6_forecast
    "src/scripts/make_paper_figures.py",          # fig_F5_calibration (+ fig_F1, fig_F3)
    "src/experiments/run_cmb_earlyde.py",         # fig_cmb_earlyde        (Fig 5)
    "src/experiments/run_e1_figure.py",           # fig_e1_mechanism       (Fig 9)
    "src/experiments/run_delta_orthogonality.py", # fig_delta_orthogonality (Fig 10)
    "src/experiments/run_chr_experiments.py",     # fig_chr_experiments
    "src/experiments/run_syk_scrambling.py",      # fig_syk_scrambling     (Fig G1)
]

for sc in SCRIPTS:
    print(f"=== {sc} ===")
    try:
        runpy.run_path(os.path.join(ROOT, sc), run_name="__main__")
    except SystemExit:
        pass
    except Exception as e:  # one bad script should not block the rest
        print(f"  {sc} FAILED: {e}")

print("\nDone. Vector PDFs written next to the PNGs in results/ and results/paper/.")
