"""Minimal serializer for OpenFOAM ASCII dictionaries.

OpenFOAM case files are C++-style dictionaries: a ``FoamFile`` header followed by
``key value;`` entries, nested ``key { ... }`` sub-dicts and ``( ... )`` lists.
This module renders Python dicts/lists into that syntax so :mod:`pyflam.cfd.case`
and :mod:`pyflam.cfd.mesh` can write cases without a PyFoam dependency.

Conventions:
    * Python ``dict``  -> a sub-dictionary block ``key { ... }``.
    * ``list``/``tuple`` -> an OpenFOAM list ``( a b c )``.
    * ``bool`` -> ``true`` / ``false``; ``int``/``float`` -> their literal text.
    * ``str`` is emitted **verbatim** — pass things like ``"uniform (0 0 0)"`` or
      ``"[0 1 -1 0 0 0 0]"`` as strings.
"""

from __future__ import annotations

_HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: written by pyflam                      |
|  \\\\    /   O peration     |                                                 |
|   \\\\  /    A nd           |                                                 |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    {location}object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""

_FOOTER = "\n// ************************************************************************* //\n"


def _fmt(value) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, str):
        return value  # verbatim
    if isinstance(value, (list, tuple)):
        return "(" + " ".join(_fmt(v) for v in value) + ")"
    raise TypeError(f"cannot serialize {type(value).__name__} to OpenFOAM syntax")


def render(entries: dict, level: int = 0) -> str:
    """Render a dict of entries to OpenFOAM dictionary text (no header)."""
    pad = "    " * level
    lines: list[str] = []
    for key, value in entries.items():
        if isinstance(value, dict):
            lines.append(f"{pad}{key}")
            lines.append(f"{pad}{{")
            lines.append(render(value, level + 1))
            lines.append(f"{pad}}}")
        else:
            gap = " " * max(1, 16 - len(key))
            lines.append(f"{pad}{key}{gap}{_fmt(value)};")
    return "\n".join(lines)


def header(cls: str, obj: str, location: str | None = None) -> str:
    """Render just the ``FoamFile`` header block."""
    loc = f'location    "{location}";\n    ' if location else ""
    return _HEADER.format(cls=cls, obj=obj, location=loc)


def foam_file(cls: str, obj: str, entries: dict,
              location: str | None = None) -> str:
    """Render a full OpenFOAM file: header + entries + footer."""
    return header(cls, obj, location) + "\n" + render(entries) + "\n" + _FOOTER


def write(path: str, cls: str, obj: str, entries: dict,
          location: str | None = None) -> None:
    """Write a full OpenFOAM dictionary file to ``path``."""
    with open(path, "w") as fh:
        fh.write(foam_file(cls, obj, entries, location))
