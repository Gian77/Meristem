#!/usr/bin/env python3
"""Create dated report and annotation archives from a completed output directory."""

import argparse
import datetime as dt
import zipfile
from pathlib import Path


def files_under(root: Path, prefixes=None):
    return [
        path for path in sorted(root.rglob("*"))
        if path.is_file() and "deliverables" not in path.relative_to(root).parts
        and (prefixes is None or any(path.relative_to(root).as_posix() == p or path.relative_to(root).as_posix().startswith(p + "/") for p in prefixes))
    ]


def build_readme(title, outdir, files, date_text, purpose):
    lines = [
        f"# {title}",
        "",
        f"Created: {date_text}",
        f"Source output directory: `{outdir}`",
        "",
        purpose, "", "The input assemblies and shared databases are intentionally not included.",
        "", "## How to use this archive", "",
        "- `pipeline_info/trace.tsv` is a machine-readable Nextflow task trace: it records task names, status, duration, CPU, memory, and I/O. It is execution provenance, not a biological quality score.",
        "- `pipeline_info/report.html`, `timeline.html`, and `dag.html` summarize execution structure and resource use.",
        "- `reports/` contains the human-readable analysis, interpretation, tools guide, and portable report formats.",
        "- `gene_prediction/` contains the GFF3 annotations and translated proteins; `busco_comparison/` contains marker-completeness evidence.",
        "- `assembly_qc/`, `repeats/`, and `comparison/` provide the assembly, repeat-masking, and structural-comparison evidence needed to judge the annotations.",
        "",
        "## Complete file listing", "",
    ]
    lines.extend(f"- `{path.relative_to(outdir).as_posix()}`" for path in files)
    lines.extend(["", "All paths above are relative to the output directory.", ""])
    return "\n".join(lines)


def package(outdir: Path, date_text: str):
    deliverables = outdir / "deliverables"
    deliverables.mkdir(parents=True, exist_ok=True)
    report_files = files_under(outdir)
    annotation_files = files_under(outdir, ["assembly_qc", "repeats", "gene_prediction", "comparison", "busco_comparison", "functional", "gene_origin"])
    archives = []
    for name, title, files, purpose in [
        (f"meristem_results_{date_text}.zip", "Meristem results package", report_files,
         "This archive contains the complete published reports, quality-control outputs, execution metadata, and annotation outputs from the run."),
        (f"meristem_annotations_{date_text}.zip", "Meristem annotation package", annotation_files,
         "This archive contains the actual structural annotations and their related assembly, repeat, comparison, and BUSCO files."),
    ]:
        archive = deliverables / name
        readme = build_readme(title, outdir, files, date_text, purpose)
        # BUSCO output contains many already-compressed/small files; storing
        # them avoids spending minutes recompressing the handoff and keeps the
        # package step reliable on the login node.
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_STORED) as handle:
            handle.writestr("README.md", readme)
            for path in files:
                handle.write(path, path.relative_to(outdir).as_posix())
        archives.append(archive)
    for old in (deliverables / "RESULTS_INDEX.md", deliverables / "meristem_annotation_deliverables.zip"):
        if old.exists():
            old.unlink()
    return archives


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--date", help="DD-MM-YYYY; defaults to the current date")
    args = parser.parse_args()
    date_text = args.date or dt.date.today().strftime("%d-%m-%Y")
    for archive in package(args.outdir.resolve(), date_text):
        print(f"Created {archive}")


if __name__ == "__main__":
    main()
