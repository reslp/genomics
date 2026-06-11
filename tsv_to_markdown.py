#!/usr/bin/env python3
"""Convert a TSV file to a Markdown table.

Usage:
  python3 tsv_to_markdown.py input.tsv > table.md
  python3 tsv_to_markdown.py input.tsv -o table.md
  python3 tsv_to_markdown.py input.tsv --no-header -o table.md
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def escape_markdown_cell(value: str) -> str:
    """Escape Markdown table separators and normalize whitespace."""
    value = value.strip().replace("\\", "\\\\").replace("|", "\\|")
    return value if value else " "


def normalize_rows(rows: list[list[str]]) -> list[list[str]]:
    """Pad rows to the same width so the Markdown table is valid."""
    if not rows:
        return []
    width = max(len(row) for row in rows)
    return [row + [""] * (width - len(row)) for row in rows]


def tsv_to_markdown(tsv_path: Path, has_header: bool = True) -> str:
    with tsv_path.open(newline="", encoding="utf-8") as handle:
        rows = normalize_rows(list(csv.reader(handle, delimiter="\t")))

    if not rows:
        return ""

    if has_header:
        header = rows[0]
        body = rows[1:]
    else:
        header = [f"Column {i}" for i in range(1, len(rows[0]) + 1)]
        body = rows

    lines = [
        "| " + " | ".join(escape_markdown_cell(cell) for cell in header) + " |",
        "| " + " | ".join("---" for _ in header) + " |",
    ]

    for row in body:
        lines.append("| " + " | ".join(escape_markdown_cell(cell) for cell in row) + " |")

    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert a TSV file to a Markdown table.")
    parser.add_argument("input", type=Path, help="Input TSV file")
    parser.add_argument("-o", "--output", type=Path, help="Output Markdown file; defaults to stdout")
    parser.add_argument("--no-header", action="store_true", help="Treat all rows as data and create generic headers")
    args = parser.parse_args()

    markdown = tsv_to_markdown(args.input, has_header=not args.no_header)

    if args.output:
        args.output.write_text(markdown, encoding="utf-8")
    else:
        sys.stdout.write(markdown)


if __name__ == "__main__":
    main()
