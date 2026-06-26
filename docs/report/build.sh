#!/bin/sh
# Build the scientific report in PDF + DOCX (English and Italian) from the markdown
# sources, using pandoc. Generated files land next to this script and are git-ignored
# (regenerate locally or in CI). Run from anywhere.
#
#   sh docs/report/build.sh
#
# Needs: pandoc, and a PDF engine (tectonic or xelatex) with a Unicode-capable serif
# font. The script auto-detects both; override with PDF_ENGINE and MAINFONT.
set -e

HERE=$(cd "$(dirname "$0")" && pwd)
DOCS=$(cd "$HERE/.." && pwd)

# --- pick a PDF engine ---
if [ -z "$PDF_ENGINE" ]; then
    if command -v tectonic >/dev/null 2>&1; then PDF_ENGINE=tectonic
    elif command -v xelatex >/dev/null 2>&1; then PDF_ENGINE=xelatex
    elif command -v lualatex >/dev/null 2>&1; then PDF_ENGINE=lualatex
    else echo "no PDF engine (tectonic / xelatex / lualatex) found" >&2; exit 1
    fi
fi

# --- pick a Unicode serif font present on this machine ---
if [ -z "$MAINFONT" ]; then
    for f in "STIX Two Text" "DejaVu Serif" "Noto Serif" "TeX Gyre Termes" "FreeSerif"; do
        if fc-list 2>/dev/null | grep -qi "$f"; then MAINFONT="$f"; break; fi
    done
    [ -z "$MAINFONT" ] && MAINFONT="DejaVu Serif"
fi
# a mono font for code blocks (box-drawing in equations)
if [ -z "$MONOFONT" ]; then
    for f in "Menlo" "DejaVu Sans Mono" "Noto Sans Mono" "FreeMono"; do
        if fc-list 2>/dev/null | grep -qi "$f"; then MONOFONT="$f"; break; fi
    done
    [ -z "$MONOFONT" ] && MONOFONT="DejaVu Sans Mono"
fi

echo "PDF engine: $PDF_ENGINE | mainfont: $MAINFONT | monofont: $MONOFONT"

build() {  # src  out_base  title
    src="$1"; out="$HERE/$2"; title="$3"
    pandoc "$src" -o "${out}.docx" --toc --toc-depth=2 -M title="$title"
    # PDF: sanitize the few LaTeX-text-mode-unfriendly glyphs first (markdown keeps them).
    python3 "$HERE/_sanitize_for_pdf.py" < "$src" \
      | pandoc -o "${out}.pdf" --pdf-engine="$PDF_ENGINE" --toc --toc-depth=2 \
          -V geometry:margin=2.2cm -V fontsize=10pt -V colorlinks=true \
          -V mainfont="$MAINFONT" -V monofont="$MONOFONT" -M title="$title"
    echo "  wrote ${out}.pdf and ${out}.docx"
}

build "$DOCS/pyflam_scientific_report.md"    pyflam_report_EN \
      "pyflam — Scientific, Technical and Operational Report"
build "$DOCS/pyflam_scientific_report_IT.md" pyflam_report_IT \
      "pyflam — Relazione Scientifica, Tecnica e Operativa"

echo "done -> $HERE/{pyflam_report_EN,pyflam_report_IT}.{pdf,docx}"
