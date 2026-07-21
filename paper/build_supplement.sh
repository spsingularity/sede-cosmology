#!/usr/bin/env bash
# Build the JCAP-class Supplementary Material PDF for Paper I from SEDE_supplement.md.
#   makedoc (title/abstract -> YAML) -> strip the citeproc refs placeholder
#   -> pandoc --natbib (emits \citep) with template_cosmology.tex -> pdflatex + bibtex.
set -e
cd "$(dirname "$0")"
mkdir -p tex

python3 tools/makedoc.py SEDE_supplement.md .SEDE.supplement.md
trap 'rm -f .SEDE.supplement.md' EXIT
# drop the citeproc "## References / ::: {#refs}" placeholder — natbib \bibliography handles it
python3 - <<'PY'
import re
t=open('.SEDE.supplement.md',encoding='utf-8').read()
t=re.sub(r'\n##\s+References\s*\n+:::\s*\{#refs\}\s*\n:::\s*\n','\n',t)
open('.SEDE.supplement.md','w',encoding='utf-8').write(t)
PY

pandoc -f markdown-superscript-subscript .SEDE.supplement.md -o tex/SEDE_supplement.tex \
  --standalone \
  --shift-heading-level-by=-1 \
  --natbib \
  --template=tools/template_cosmology.tex

# figure paths: from tex/ the repo results dir is ../../results ; prefer vector .pdf
perl -0pi -e 's#\{\.\./results/([^}]+?)\.png\}#{../../results/\1}#g' tex/SEDE_supplement.tex
# numeric journal (JHEP): no author-prominent form — normalise \citet -> \citep (silences natbib)
perl -0pi -e 's#\\citet\{#\\citep{#g' tex/SEDE_supplement.tex

( cd tex && \
  pdflatex -interaction=nonstopmode SEDE_supplement.tex >SEDE_supplement.build.log 2>&1 ; \
  BIBINPUTS="..:$BIBINPUTS" bibtex SEDE_supplement      >>SEDE_supplement.build.log 2>&1 ; \
  pdflatex -interaction=nonstopmode SEDE_supplement.tex >>SEDE_supplement.build.log 2>&1 ; \
  pdflatex -interaction=nonstopmode SEDE_supplement.tex >>SEDE_supplement.build.log 2>&1 ) || true

if [ -f tex/SEDE_supplement.pdf ]; then
  echo "built paper/tex/SEDE_supplement.pdf"
  cp tex/SEDE_supplement.pdf SEDE_supplement.pdf     # keep the final PDF alongside the manuscript
  grep -c "^!" tex/SEDE_supplement.build.log | awk '{print $1" LaTeX errors (see tex/SEDE_supplement.build.log)"}'
  grep -c "Warning--" tex/SEDE_supplement.build.log 2>/dev/null | awk '{print $1" bibtex warnings"}'
else
  echo "BUILD FAILED — see tex/SEDE_supplement.build.log"; grep -A2 '^!' tex/SEDE_supplement.build.log | head -20
fi
