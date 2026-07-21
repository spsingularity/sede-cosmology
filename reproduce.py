#!/usr/bin/env python3
"""Single entry point that regenerates the paper's figures from the vendored code.

    python reproduce.py

What runs here (self-contained: the vendored `sede` package + the small data vectors in
`data/`):
  figures/figures.py           -> Figs 1-8 : fig1_mechanism, fig2_inputs, fig3_datafit,
                                             fig4_eos, fig5_oos, fig6_forecast  (results/paper/)
  src/experiments/run_cmb_earlyde.py        -> Fig 5  : fig_cmb_earlyde        (results/)
  src/experiments/run_delta_orthogonality.py-> Fig 10 : fig_delta_orthogonality(results/)

Two figures need heavier inputs that live with the code repository (see data/README.md),
so they are shipped pre-built and only regenerated when those inputs are present:
  src/experiments/run_e1_figure.py -> Fig 9  : fig_e1_mechanism    (needs CAMB: pip install camb)
  fig_F5_calibration           -> Fig 6  : false-preference null   (needs the mock-calibration
                                           driver + the full Pantheon+ covariance)

Notes:
  * fig3_datafit falls back to the stat-only Pantheon+ errors when the 32 MB
    Pantheon+SH0ES_STAT+SYS.cov is absent (a WARN is printed); drop that file into data/ for
    the full covariance (data/README.md).
"""
import os, sys, subprocess

ROOT = os.path.dirname(os.path.abspath(__file__))
_paths = [os.path.join(ROOT, "src"), os.path.join(ROOT, "src", "experiments")]
if os.environ.get("PYTHONPATH"):
    _paths.append(os.environ["PYTHONPATH"])
ENV = dict(os.environ, PYTHONPATH=os.pathsep.join(_paths))


def run(cmd, optional=False):
    print(f"\n=== {' '.join(cmd)} ===", flush=True)
    r = subprocess.run([sys.executable, *cmd], cwd=ROOT, env=ENV)
    if r.returncode and not optional:
        print(f"  FAILED ({r.returncode})")
    return r.returncode == 0


def have_camb():
    try:
        import camb  # noqa: F401
        return True
    except Exception:
        return False


def main():
    run(["figures/figures.py"])
    run(["src/experiments/run_delta_orthogonality.py"])   # imported by run_cmb_earlyde; run first
    run(["src/experiments/run_cmb_earlyde.py"])
    if have_camb():
        run(["src/experiments/run_e1_figure.py"], optional=True)
    else:
        print("\n[skip] fig_e1_mechanism needs CAMB (pip install camb) — shipped pre-built.")
    print("\n[note] fig_F5_calibration is the mock false-preference run; it needs the "
          "calibration driver and the full Pantheon+ covariance (data/README.md) and is "
          "shipped pre-built.")
    print("\nFigures are in results/ and results/paper/.")


if __name__ == "__main__":
    main()
