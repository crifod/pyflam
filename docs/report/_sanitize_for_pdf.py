#!/usr/bin/env python
"""Map the few LaTeX-text-mode-unfriendly Unicode symbols in a markdown file to
readable ASCII/word equivalents, for PDF typesetting via pandoc + tectonic.

Greek letters and sub/superscripts render fine in STIX Two Text, so they are kept;
only arrows, relations and operators that the xelatex tex-text mapping intercepts
are transliterated. Reads stdin, writes stdout. The source markdown is unchanged.
"""
import sys

REPL = {
    "→": " -> ", "⇒": " => ", "⇄": " <-> ",
    "∝": " proportional to ", "∈": " in ",
    "≈": " ~ ", "≠": " != ", "≤": " <= ", "≥": " >= ",
    "∇": "grad ", "∂": "d", "∑": "sum", "√": "sqrt",
    "×": "x", "÷": "/", "±": "+/-", "½": "1/2",
    "′": "'", "·": "*", "─": "-",
    "̂": "",            # combining circumflex (k-hat) -> drop
    "k̂": "k_hat",
}


def main():
    text = sys.stdin.read()
    text = text.replace("k̂", "k_hat")
    for a, b in REPL.items():
        text = text.replace(a, b)
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
