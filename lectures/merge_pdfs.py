#!/usr/bin/env python3
"""
Concatenate every PDF found under a directory (recursively) into one PDF.

Usage:
    python merge_pdfs.py ROOT_DIR [-o merged.pdf] [--no-bookmarks] [--dry-run]

Files are ordered by folder, then by filename, using "natural" sorting
(so "file2.pdf" comes before "file10.pdf"). Each source file gets a
bookmark in the output, nested under its folder, so you can navigate easily.

Requires: pip install pypdf
"""

import argparse
import re
import sys
from pathlib import Path

from pypdf import PdfReader, PdfWriter


def natural_key(path: Path, root: Path):
    """Sort key that treats digit runs as numbers, compared part by part."""
    parts = path.relative_to(root).parts
    return [
        [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", part)]
        for part in parts
    ]


def find_pdfs(root: Path, exclude: Path):
    pdfs = [
        p for p in root.rglob("*")
        if p.is_file() and p.suffix.lower() == ".pdf" and p.resolve() != exclude
    ]
    return sorted(pdfs, key=lambda p: natural_key(p, root))


def main():
    parser = argparse.ArgumentParser(description="Merge all PDFs in subdirectories into one PDF.")
    parser.add_argument("root", type=Path, help="Top-level directory to search")
    parser.add_argument("-o", "--output", type=Path, default=Path("merged.pdf"),
                        help="Output file (default: merged.pdf)")
    parser.add_argument("--no-bookmarks", action="store_true",
                        help="Don't add a bookmark for each source file")
    parser.add_argument("--dry-run", action="store_true",
                        help="Just list the files in the order they'd be merged")
    args = parser.parse_args()

    root = args.root.resolve()
    if not root.is_dir():
        sys.exit(f"Error: {root} is not a directory")

    output = args.output.resolve()
    pdfs = find_pdfs(root, exclude=output)
    if not pdfs:
        sys.exit(f"No PDF files found under {root}")

    if args.dry_run:
        for p in pdfs:
            print(p.relative_to(root))
        print(f"\n{len(pdfs)} file(s) would be merged into {output}")
        return

    writer = PdfWriter()
    folder_marks = {}   # folder path -> bookmark object
    merged, skipped = 0, []

    for pdf in pdfs:
        rel = pdf.relative_to(root)
        try:
            reader = PdfReader(pdf)
            if reader.is_encrypted:
                reader.decrypt("")  # works for PDFs with an empty user password
            start_page = len(writer.pages)
            for page in reader.pages:
                writer.add_page(page)
        except Exception as e:
            skipped.append((rel, e))
            print(f"  ! Skipped {rel}: {e}", file=sys.stderr)
            continue

        if not args.no_bookmarks and len(writer.pages) > start_page:
            parent = None
            folder = rel.parent
            if folder != Path("."):
                if folder not in folder_marks:
                    folder_marks[folder] = writer.add_outline_item(str(folder), start_page)
                parent = folder_marks[folder]
            writer.add_outline_item(pdf.stem, start_page, parent=parent)

        merged += 1
        print(f"  + {rel} ({len(reader.pages)} pages)")

    if merged == 0:
        sys.exit("Nothing could be merged.")

    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "wb") as f:
        writer.write(f)

    print(f"\nMerged {merged} file(s), {len(writer.pages)} pages total -> {output}")
    if skipped:
        print(f"Skipped {len(skipped)} file(s) that couldn't be read.")


if __name__ == "__main__":
    main()
