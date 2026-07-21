#!/usr/bin/env bash
# Build the JCAP-class PDF (jcappub.sty) for Paper I from SEDE_cosmology.md.
#   makedoc (title/abstract -> YAML) -> strip the citeproc refs placeholder
#   -> pandoc --natbib (emits \citep) with template_cosmology.tex -> xelatex + bibtex.
set -e
cd "$(dirname "$0")"
mkdir -p tex

python3 tools/makedoc.py SEDE_cosmology.md .SEDE.cosmology.md
trap 'rm -f .SEDE.cosmology.md' EXIT
# drop the citeproc "## References / ::: {#refs}" placeholder — natbib \bibliography handles it
python3 - <<'PY'
import re
t=open('.SEDE.cosmology.md',encoding='utf-8').read()
t=re.sub(r'\n##\s+References\s*\n+:::\s*\{#refs\}\s*\n:::\s*\n','\n',t)
open('.SEDE.cosmology.md','w',encoding='utf-8').write(t)
PY

pandoc -f markdown-superscript-subscript .SEDE.cosmology.md -o tex/SEDE_cosmology.tex \
  --standalone \
  --shift-heading-level-by=-1 \
  --natbib \
  --template=tools/template_cosmology.tex

# figure paths: from tex/ the repo results dir is ../../results ; prefer vector .pdf
perl -0pi -e 's#\{\.\./results/([^}]+?)\.png\}#{../../results/\1}#g' tex/SEDE_cosmology.tex
# numeric journal (JHEP): no author-prominent form — normalise \citet -> \citep (silences natbib)
perl -0pi -e 's#\\citet\{#\\citep{#g' tex/SEDE_cosmology.tex
# number display equations: pandoc emits \[ ... \]; convert to a numbered equation environment
perl -0pi -e 's/\\\[(.*?)\\\]/\\begin{equation}$1\\end{equation}/gs' tex/SEDE_cosmology.tex

( cd tex && \
  pdflatex -interaction=nonstopmode SEDE_cosmology.tex >SEDE_cosmology.build.log 2>&1 ; \
  BIBINPUTS="..:$BIBINPUTS" bibtex SEDE_cosmology      >>SEDE_cosmology.build.log 2>&1 ; \
  pdflatex -interaction=nonstopmode SEDE_cosmology.tex >>SEDE_cosmology.build.log 2>&1 ; \
  pdflatex -interaction=nonstopmode SEDE_cosmology.tex >>SEDE_cosmology.build.log 2>&1 ) || true

if [ -f tex/SEDE_cosmology.pdf ]; then
  echo "built paper/tex/SEDE_cosmology.pdf"
  cp tex/SEDE_cosmology.pdf SEDE_cosmology.pdf       # keep the committed final PDF in sync
  grep -c "^!" tex/SEDE_cosmology.build.log | awk '{print $1" LaTeX errors (see tex/SEDE_cosmology.build.log)"}'
  grep -c "Warning--" tex/SEDE_cosmology.build.log 2>/dev/null | awk '{print $1" bibtex warnings"}'
else
  echo "BUILD FAILED — see tex/SEDE_cosmology.build.log"; grep -A2 '^!' tex/SEDE_cosmology.build.log | head -20
fi
