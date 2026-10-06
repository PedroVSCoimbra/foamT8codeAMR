#!/usr/bin/env python3
"""Doxygen input filter for OpenFOAM-style sources (.H/.C).

Replaces foam-extend's bin/doxyFilt, whose License pattern targets the old
GPL text and therefore drops the whole Description of GPLv3 headers.

Transformations (the line count is preserved, so source browsing and
warning line numbers match the original file):

- The header banner (License, Class, Description, Usage, SourceFiles, ...)
  becomes a Doxygen block: Class -> \\class, the first paragraph of
  Description -> \\brief and the rest -> detailed description, Usage and
  SourceFiles -> paragraphs. The License text is dropped.
- "//-" comments, and the "//  " lines continuing them, become "//!" so
  that Doxygen attaches them to the following declaration.

Usage (as Doxygen INPUT_FILTER / FILTER_PATTERNS):
    python3 doxyFilt.py <file>
"""

import os
import re
import sys

SECTION = re.compile(
    r"^(License|Class|Description|Usage|SourceFiles|See also|SeeAlso|"
    r"Author|Note|Application|InClass|Typedef|Global)\s*$"
)
BANNER_START = re.compile(r"^/\*-+\*\\\s*$")
BANNER_END = re.compile(r"^\\\*-+\*/\s*$")
DOC_START = re.compile(r"^(\s*)//-\s?(.*)$")
DOC_CONT = re.compile(r"^(\s*)//  ?(.*)$")


def dedent(lines, n=4):
    return [ln[n:] if ln.startswith(" " * n) else ln.lstrip() for ln in lines]


def split_sections(body):
    """Ordered {section name: lines} of a banner body."""
    sections, current = {}, None
    for line in body:
        match = SECTION.match(line)
        if match:
            current = match.group(1)
            sections[current] = []
        elif current is not None:
            sections[current].append(line)
    return {k: dedent(v) for k, v in sections.items()}


def strip_blank(lines):
    while lines and not lines[0].strip():
        lines = lines[1:]
    while lines and not lines[-1].strip():
        lines = lines[:-1]
    return lines


def banner_to_doxygen(body, file_name):
    """Doxygen comment lines replacing the banner body."""
    sections = split_sections(body)
    out = [f"/*! \\file {file_name}"]

    if "Class" in sections:
        names = [ln.strip() for ln in sections["Class"] if ln.strip()]
        if names:
            out.append(f"\\class {names[0]}")

    desc = strip_blank(sections.get("Description", []))
    if desc:
        try:
            first_blank = next(i for i, ln in enumerate(desc) if not ln.strip())
        except StopIteration:
            first_blank = len(desc)
        out.append("\\brief " + " ".join(ln.strip()
                                          for ln in desc[:first_blank]))
        out.append("")
        out.extend(desc[first_blank + 1:])

    for name, title in (("Usage", "Usage"), ("Note", "Note"),
                        ("See also", "See also"), ("SeeAlso", "See also")):
        text = strip_blank(sections.get(name, []))
        if text:
            out += ["", f"\\par {title}", ""] + text

    files = [ln.strip() for ln in sections.get("SourceFiles", [])
             if ln.strip()]
    if files:
        out += ["", "\\par Source files", ""] + [f"- {f}" for f in files]

    out.append("*/")
    return out


def filter_lines(lines, file_name):
    out = []
    i = 0
    # Header banner: keep the logo lines before it, convert its body
    if lines and BANNER_START.match(lines[0]):
        end = next((k for k, ln in enumerate(lines) if k > 0
                    and BANNER_END.match(ln)), None)
        if end is not None:
            start = next((k for k in range(1, end)
                          if SECTION.match(lines[k])), end)
            # Logo part stays a plain comment, closed before the Doxygen block
            out.extend(lines[:start - 1])
            out.append(lines[start - 1].rstrip("\n") + " */\n"
                       if start > 1 else "\n")
            doxy = banner_to_doxygen(
                [ln.rstrip("\n") for ln in lines[start:end]], file_name)
            n_orig = end - start + 1
            # Preserve the line count with blank lines before the closing */
            # (a longer block only shifts the line numbers that follow)
            pad = max(0, n_orig - len(doxy))
            doxy = doxy[:-1] + [""] * pad + doxy[-1:]
            out.extend(ln + "\n" for ln in doxy)
            i = end + 1

    in_doc = False
    for line in lines[i:]:
        start = DOC_START.match(line)
        if start:
            out.append(f"{start.group(1)}//! {start.group(2)}\n")
            in_doc = True
            continue
        cont = DOC_CONT.match(line) if in_doc else None
        if cont:
            out.append(f"{cont.group(1)}//! {cont.group(2)}\n")
            continue
        in_doc = False
        out.append(line)
    return out


def main():
    path = sys.argv[1]
    with open(path, errors="replace") as fh:
        lines = fh.readlines()
    sys.stdout.writelines(filter_lines(lines, os.path.basename(path)))


if __name__ == "__main__":
    main()
