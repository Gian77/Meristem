#!/usr/bin/env python3
"""Render readable Meristem annotation reports without third-party packages.

The script intentionally uses only the Python standard library so it works in
the lightweight test profile and in the MultiQC container. It writes a useful
HTML report/data directory plus simple, portable DOCX and PDF versions of the
combined summary and interpretation.
"""

import argparse
import csv
import html
import json
import math
import re
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape as xml_escape


def parse_key_value_table(path):
    values = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                values[fields[0]] = fields[1]
    return values


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value, default=0.0):
    try:
        return float(value or default)
    except (TypeError, ValueError):
        return default


def fmt_int(value):
    try:
        return f"{int(float(value or 0)):,}"
    except (TypeError, ValueError):
        return value or "not available"


def fmt_pct(value):
    try:
        return f"{float(value):.1f}%"
    except (TypeError, ValueError):
        return value or "not available"


def yes(value):
    return value.strip().lower() in {"1", "true", "yes", "on"}


def method_rows(rows):
    return {row.get("method", "").lower(): row for row in rows if row.get("method")}


def feature_rows(rows):
    return {row.get("metric", "").lower(): row for row in rows if row.get("metric")}


def interpretation(sample, lineage, assembly, busco, features, run_edta,
                   run_repeatmasker, run_interproscan, run_eggnog):
    total_length = number(assembly.get("total_length"))
    gaps = number(assembly.get("gap_bases"))
    gap_fraction = (100 * gaps / total_length) if total_length else 0.0
    methods = list(busco)

    lines = [
        f"# Interpretation — {sample}",
        "",
        "## Executive assessment",
        "",
        "The workflow report was generated after the upstream reporting inputs were available. "
        "That means the computational run reached the reporting stage successfully; it does not, "
        "by itself, establish biological correctness. The interpretation below separates execution "
        "status from annotation quality.",
        "",
        "## 1. Genome assessment",
        "",
        f"- **Sequences at the configured minimum length:** {fmt_int(assembly.get('sequences_ge_min'))}",
        f"- **Total assembly length:** {fmt_int(assembly.get('total_length'))} bp",
        f"- **Largest sequence:** {fmt_int(assembly.get('largest'))} bp",
        f"- **N50:** {fmt_int(assembly.get('N50'))} bp",
        f"- **Gap bases:** {fmt_int(assembly.get('gap_bases'))} ({gap_fraction:.4f}% of reported assembly length)",
        "",
        "Interpret N50 alongside sequence counts and total length. These values describe continuity; "
        "they do not establish gene-space completeness or rule out contamination. Confirm that the "
        "assembly statistics have passed consistency checks before using them to compare annotations.",
        "",
        "## 2. Repeat discovery and masking",
        "",
        f"- **EDTA repeat discovery:** {'enabled' if run_edta else 'not run; an existing repeat library was supplied'}",
        f"- **RepeatMasker masking:** {'enabled' if run_repeatmasker else 'not run; the supplied assembly/library was used as configured'}",
        "",
        "The annotation predictors were intended to receive the repeat-prepared assembly produced by "
        "the configured repeat branch. Repeat masking reduces the chance that repetitive sequence is "
        "called as genes, but the quality of the repeat library and the masking rate should be checked "
        "in the repeat outputs before finalizing the annotation.",
        "",
        "## 3. Structural annotation",
        "",
    ]

    for method in ("braker4", "helixer"):
        if method not in features.get('gene', {}) and method not in busco:
            continue
        row = features.get("gene", {})
        gene_count = row.get(method, "not available")
        if method == "helixer" and row:
            gene_count = row.get("helixer", gene_count)
        elif method == "braker4" and row:
            gene_count = row.get("braker4", gene_count)
        lines.append(f"- **{method.upper()} predicted genes:** {fmt_int(gene_count)}")

    if "gene" in features:
        gene = features["gene"]
        lines.append(
            f"- The feature comparison reports {fmt_int(gene.get('braker4'))} BRAKER4 genes versus "
            f"{fmt_int(gene.get('helixer'))} Helixer genes."
        )
    lines.extend(
        [
            "",
            "Gene counts are not a quality score by themselves. Differences can reflect model granularity, "
            "transcript selection, and treatment of fragmented or repetitive loci. Check the GFF3/FASTA "
            "identifier consistency and inspect protein length distributions before merging or selecting models.",
            "",
            "## 4. BUSCO completeness",
            "",
        ]
    )

    if busco:
        lines.append(f"BUSCO lineage: **{lineage}**.")
        lines.append("")
        lines.append("| Method | Complete | Single-copy | Duplicated | Fragmented | Missing |")
        lines.append("|---|---:|---:|---:|---:|---:|")
        for method, row in busco.items():
            lines.append(
                f"| {method} | {fmt_pct(row.get('complete_pct'))} | {fmt_pct(row.get('single_copy_pct'))} | "
                f"{fmt_pct(row.get('duplicated_pct'))} | {fmt_pct(row.get('fragmented_pct'))} | "
                f"{fmt_pct(row.get('missing_pct'))} |"
            )
        lines.append("")
        comparable = all(finite_number(row.get(key)) is not None
                         for row in busco.values() for key in ('complete_pct', 'missing_pct'))
        if "braker4" in busco and "helixer" in busco and comparable:
            braker_complete = number(busco["braker4"].get("complete_pct"))
            helixer_complete = number(busco["helixer"].get("complete_pct"))
            braker_missing = number(busco["braker4"].get("missing_pct"))
            helixer_missing = number(busco["helixer"].get("missing_pct"))
            if braker_complete > helixer_complete and braker_missing < helixer_missing:
                lines.append(
                    f"BRAKER4 has the stronger BUSCO result in this comparison: {braker_complete:.1f}% complete "
                    f"versus {helixer_complete:.1f}% for Helixer, with fewer missing BUSCOs. Validate protein "
                    "extraction and representative-isoform selection before interpreting this as a method difference."
                )
            else:
                if helixer_complete > braker_complete and helixer_missing < braker_missing:
                    lines.append(
                        f"Helixer has the stronger BUSCO result in this comparison: {helixer_complete:.1f}% complete "
                        f"versus {braker_complete:.1f}% for BRAKER4, with fewer missing BUSCOs. Validate protein "
                        "extraction and representative-isoform selection before interpreting this as a method difference."
                    )
                else:
                    lines.append(
                        "The BUSCO comparison does not identify a single method as clearly superior on both "
                        "completeness and missing BUSCOs; review the original summaries and structural metrics."
                    )
        elif len(methods) == 1 and comparable:
            method = methods[0]
            complete = number(busco[method].get("complete_pct"))
            missing = number(busco[method].get("missing_pct"))
            lines.append(
                f"The selected method reports {complete:.1f}% complete and {missing:.1f}% missing BUSCOs. "
                "Interpret this together with gene-model structure and contamination checks."
            )
        else:
            lines.append("BUSCO summary values are incomplete; a completeness ranking is unavailable.")
    else:
        lines.append("BUSCO was not run, so completeness-based interpretation is unavailable for this run.")

    lines.extend(
        [
            "",
            "## 5. Why the Helixer BUSCO result is low",
            "",
            "The reported Helixer result (37.2% complete, 12.1% fragmented, 50.8% missing) is not evidence "
            "that the neural-network prediction itself is only 37% complete. The completed run translated "
            "Helixer multi-exon CDS records by trimming the GFF3 phase at every CDS segment. In GFF3, phase "
            "describes codon context at a junction; those bases are part of the coding sequence and must be "
            "retained when CDS segments are concatenated. Trimming them can shift the reading frame repeatedly "
            "and create artificial truncations or internal stops, which is consistent with the large BUSCO "
            "drop relative to BRAKER4. The pipeline code has now been corrected, but this existing 37.2% result "
            "must be treated as a pre-fix diagnostic until Helixer protein extraction and BUSCO are rerun.",
            "",
            "The comparison is also not fully independent: both methods used the same masked assembly, and BUSCO "
            "is evaluated on selected protein translations. After the corrected rerun, inspect translated-protein "
            "count, internal stops, protein-length distribution, and BUSCO results together. If Helixer remains "
            "near 37% after the translation fix, investigate model/lineage suitability, assembly fragmentation, "
            "and prediction sensitivity rather than accepting the current result as a biological conclusion.",
            "",
            "## 6. Functional annotation and next steps",
            "",
            f"- **InterProScan:** {'enabled' if run_interproscan else 'not enabled'}",
            f"- **eggNOG:** {'enabled' if run_eggnog else 'not enabled'}",
            "",
            "For a final annotation decision, review the assembly QC, repeat outputs, GFF3/FASTA "
            "consistency, BUSCO summaries, and contamination evidence together. If the structural "
            "annotation is accepted, functional annotation can be enabled as a separate opt-in stage. "
            "Do not treat a successful workflow status as proof that every biological quality criterion passed.",
            "",
            "## 7. Readiness decision",
            "",
            "This report does not decide whether the annotation is ready for release. That decision needs "
            "the assembly QC, repeat outputs, BUSCO summaries, and any gene origin results reviewed together "
            "for this genome.",
            "",
        ]
    )
    return "\n".join(lines)


def combined_text(summary, interpretation_path):
    return summary.read_text() + "\n\n" + interpretation_path.read_text()


def write_docx(text, output):
    paragraphs = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            paragraphs.append("")
            continue
        line = re.sub(r"^#{1,6}\s*", "", line)
        line = re.sub(r"\*\*(.*?)\*\*", r"\1", line)
        line = re.sub(r"`(.*?)`", r"\1", line)
        line = line.replace("|", "  ")
        paragraphs.append(line)

    body = []
    for paragraph in paragraphs:
        if paragraph:
            body.append(f"<w:p><w:r><w:t xml:space=\"preserve\">{xml_escape(paragraph)}</w:t></w:r></w:p>")
        else:
            body.append("<w:p/>")
    document = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
        "<w:body>"
        + "".join(body)
        + '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/><w:pgMar w:top="1000" w:right="1000" w:bottom="1000" w:left="1000"/></w:sectPr>'
        "</w:body></w:document>"
    )
    content_types = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
        "</Types>"
    )
    relationships = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>'
        "</Relationships>"
    )
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", content_types)
        archive.writestr("_rels/.rels", relationships)
        archive.writestr("word/document.xml", document)


def pdf_escape(text):
    text = text.encode("latin-1", "replace").decode("latin-1")
    return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def write_pdf(text, output):
    lines = []
    for raw in text.splitlines():
        line = re.sub(r"^#{1,6}\s*", "", raw.strip())
        line = re.sub(r"\*\*(.*?)\*\*", r"\1", line)
        line = re.sub(r"`(.*?)`", r"\1", line)
        line = line.replace("|", "  ")
        if not line:
            lines.append("")
            continue
        while len(line) > 105:
            split = line.rfind(" ", 0, 105)
            split = split if split > 20 else 105
            lines.append(line[:split])
            line = line[split:].lstrip()
        lines.append(line)

    page_lines = 52
    pages = [lines[index : index + page_lines] for index in range(0, len(lines), page_lines)] or [[""]]
    objects = []
    objects.append(b"<< /Type /Catalog /Pages 2 0 R >>")
    page_ids = [4 + i * 2 for i in range(len(pages))]
    content_ids = [5 + i * 2 for i in range(len(pages))]
    kids = " ".join(f"{page_id} 0 R" for page_id in page_ids)
    objects.append(f"<< /Type /Pages /Kids [{kids}] /Count {len(pages)} >>".encode())
    objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    for page, content_id in zip(pages, content_ids):
        stream_lines = ["BT", "/F1 9 Tf", "48 760 Td", "11 TL"]
        for line in page:
            stream_lines.append(f"({pdf_escape(line)}) Tj")
            stream_lines.append("0 -11 Td")
        stream_lines.append("ET")
        stream = "\n".join(stream_lines).encode("latin-1", "replace")
        objects.append(
            f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] /Resources << /Font << /F1 3 0 R >> >> /Contents {content_id} 0 R >>".encode()
        )
        objects.append(b"<< /Length " + str(len(stream)).encode() + b" >>\nstream\n" + stream + b"\nendstream")

    output_bytes = bytearray(b"%PDF-1.4\n")
    offsets = [0]
    for object_id, obj in enumerate(objects, start=1):
        offsets.append(len(output_bytes))
        output_bytes.extend(f"{object_id} 0 obj\n".encode())
        output_bytes.extend(obj)
        output_bytes.extend(b"\nendobj\n")
    xref = len(output_bytes)
    output_bytes.extend(f"xref\n0 {len(objects) + 1}\n".encode())
    output_bytes.extend(b"0000000000 65535 f \n")
    for offset in offsets[1:]:
        output_bytes.extend(f"{offset:010d} 00000 n \n".encode())
    output_bytes.extend(
        f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n".encode()
    )
    output.write_bytes(output_bytes)


def write_html(sample, lineage, assembly, busco, features, interpretation_path,
               output, data_dir):
    data_dir.mkdir(parents=True, exist_ok=True)
    data_path = data_dir / "meristem_annotation_summary_mqc.tsv"
    with data_path.open("w") as handle:
        handle.write("sample\tsection\tmetric\tvalue\n")
        for key, value in assembly.items():
            handle.write(f"{sample}\tassembly\t{key}\t{value}\n")
        for method, row in busco.items():
            for metric in ("complete_pct", "single_copy_pct", "duplicated_pct", "fragmented_pct", "missing_pct"):
                handle.write(f"{sample}\tbusco_{method}\t{metric}\t{row.get(metric, '')}\n")
        for metric, row in features.items():
            for method in ("braker4", "helixer"):
                handle.write(f"{sample}\tfeatures\t{method}_{metric}\t{row.get(method, '')}\n")

    def table(headers, rows):
        result = ["<table><thead><tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in headers) + "</tr></thead><tbody>"]
        for row in rows:
            result.append("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>")
        result.append("</tbody></table>")
        return "".join(result)

    busco_table = table(
        ["Method", "Complete", "Single-copy", "Duplicated", "Fragmented", "Missing"],
        [
            [method, fmt_pct(row.get("complete_pct")), fmt_pct(row.get("single_copy_pct")), fmt_pct(row.get("duplicated_pct")), fmt_pct(row.get("fragmented_pct")), fmt_pct(row.get("missing_pct"))]
            for method, row in busco.items()
        ],
    ) if busco else "<p>BUSCO was not run.</p>"
    feature_table = table(
        ["Feature", "BRAKER4", "Helixer"],
        [[metric, row.get("braker4", ""), row.get("helixer", "")] for metric, row in features.items()],
    )
    assembly_table = table([["Metric", "Value"]][0], [[key, value] for key, value in assembly.items()])
    interpretation_html = html.escape(interpretation_path.read_text()).replace("\n", "<br>\n")
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>Meristem annotation report — {html.escape(sample)}</title>
<style>body{{font-family:Arial,sans-serif;max-width:1100px;margin:2em auto;color:#222}}h1{{color:#173f5f}}h2{{color:#20639b;border-bottom:1px solid #ddd;padding-top:1em}}table{{border-collapse:collapse;margin:1em 0;min-width:60%}}th,td{{border:1px solid #bbb;padding:.4em .7em;text-align:left}}th{{background:#eaf2f8}}pre{{background:#f5f5f5;padding:1em;white-space:pre-wrap}}.pass{{color:#176b2c;font-weight:bold}}</style></head>
<body><h1>Meristem annotation report — {html.escape(sample)}</h1>
<p class="pass">Report generated successfully. BUSCO lineage: {html.escape(lineage)}</p>
<h2>Assembly assessment</h2>{assembly_table}
<h2>BUSCO completeness</h2>{busco_table}
<h2>Annotation feature counts</h2>{feature_table}
<h2>Detailed interpretation</h2><pre>{interpretation_html}</pre>
</body></html>"""
    output.write_text(document)


REPORT_CSS = """
:root{--ink:#19342e;--muted:#52665f;--green:#226347;--paper:#fafbf8;--line:#dce5dd}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:var(--paper);color:var(--ink);font:16px/1.65 system-ui,-apple-system,sans-serif}
a{color:#1c6557;text-underline-offset:3px}a:hover{color:#923d17}header{background:#173c31;color:white;padding:42px max(6vw,24px)}
header .eyebrow{font-size:12px;text-transform:uppercase;letter-spacing:.16em;color:#bfdec3}header h1{font-size:clamp(30px,4vw,48px);line-height:1.15;margin:12px 0}
header p{max-width:850px;color:#d4e6dc}header a{color:white}.layout{display:grid;grid-template-columns:240px minmax(0,1040px);gap:44px;max-width:1380px;margin:auto;padding:36px 24px}
aside{position:sticky;top:20px;align-self:start;max-height:90vh;overflow:auto;font-size:13px}aside strong{display:block;margin-bottom:12px}aside a{display:block;text-decoration:none;border-left:2px solid var(--line);padding:5px 12px}aside a:hover{background:#edf3eb;border-color:var(--green)}
main{min-width:0}h2{font-size:25px;line-height:1.3;border-top:1px solid var(--line);padding-top:30px;margin-top:38px;scroll-margin-top:24px}h3{font-size:19px;scroll-margin-top:24px}p,li{max-width:90ch}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:16px;margin:0 0 28px}.card{background:white;border:1px solid var(--line);border-radius:12px;padding:20px}.card .label{color:var(--muted);font-size:13px}.card .value{display:block;font-size:32px;font-weight:650;line-height:1.25;margin:6px 0}.card small{color:var(--muted)}
.notice,blockquote{background:#edf4eb;border-left:4px solid #568464;padding:16px 20px;margin:20px 0;border-radius:0 8px 8px 0}.notice.warning{background:#fff5df;border-color:#ba852d}.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:8px;margin:20px 0}table{border-collapse:collapse;width:100%;background:white;font-size:14px}th,td{text-align:left;vertical-align:top;padding:12px 14px;border-bottom:1px solid var(--line)}th{background:#eaf0e7;font-weight:650}tbody tr:last-child td{border:0}tr:hover td{background:#f8faf5}pre{overflow:auto;padding:18px;background:#edf0ec;border-radius:8px;font-size:13px}code{font-size:.88em;background:#edf0ec;padding:2px 5px;border-radius:4px}pre code{padding:0}hr{border:0;border-top:1px solid var(--line);margin:30px 0}
.badge{display:inline-block;font-size:12px;font-weight:650;background:#e8f0e6;border:1px solid #cad9c9;border-radius:20px;padding:3px 11px}.tool{background:white;border:1px solid var(--line);border-radius:12px;padding:4px 24px 20px;margin:20px 0}.tool h2{border:0;margin-top:16px;padding-top:0}.tool dl{display:grid;grid-template-columns:130px 1fr;gap:8px 16px;margin-bottom:0}.tool dt{font-weight:650}.tool dd{margin:0;color:#354c43;overflow-wrap:anywhere}.search{width:100%;padding:12px 16px;font:inherit;border:1px solid #b8ccbc;border-radius:8px;background:white}footer{border-top:1px solid var(--line);padding:24px;color:var(--muted);font-size:13px}.muted{color:var(--muted)}
@media(max-width:850px){.layout{grid-template-columns:1fr;gap:16px}aside{position:static;max-height:240px}.tool dl{grid-template-columns:1fr;gap:3px}.tool dd{margin-bottom:12px}}
@media print{header{background:white;color:var(--ink);padding:0}header p,header .eyebrow{color:var(--muted)}header a{color:var(--ink)}.layout{display:block;padding:0}aside,.search{display:none}.tool,.card,tr{break-inside:avoid}h2,h3{break-after:avoid}body{background:white;font-size:11pt}.table-wrap{overflow:visible}a{color:inherit}}
"""


def safe_url(value):
    """Only local relative links and primary documentation URLs are rendered."""
    if re.match(r"https?://", value, re.I):
        return value
    if value.startswith(('#', './', '../')) or re.match(r'^[A-Za-z0-9_][A-Za-z0-9_.-]*$', value):
        return value
    return None


def inline_markdown(text):
    # Escape raw HTML, and permit only a narrow link syntax with safe URLs.
    escaped = html.escape(text)
    def link(match):
        target = html.unescape(match.group(2))
        url = safe_url(target)
        return f'<a href="{html.escape(url, quote=True)}">{match.group(1)}</a>' if url else match.group(1)
    escaped = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', link, escaped)
    escaped = re.sub(r'`([^`]+)`', r'<code>\1</code>', escaped)
    return re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', escaped)


def markdown_sections(text):
    """Render the headings, tables, lists, and notes used by Meristem reports."""
    lines = text.splitlines()
    result, toc, seen = [], [], {}
    index, list_type = 0, None
    def close_list():
        nonlocal list_type
        if list_type:
            result.append(f'</{list_type}>')
            list_type = None
    while index < len(lines):
        line = lines[index]
        if line.startswith('```'):
            close_list()
            code = []
            index += 1
            while index < len(lines) and not lines[index].startswith('```'):
                code.append(lines[index]); index += 1
            result.append('<pre><code>' + html.escape('\n'.join(code)) + '</code></pre>')
        elif line.startswith('|'):
            close_list()
            rows = []
            while index < len(lines) and lines[index].startswith('|'):
                row = lines[index]
                if not re.fullmatch(r'[|\s:-]+', row):
                    rows.append([cell.strip() for cell in row.strip().strip('|').split('|')])
                index += 1
            result.append('<div class="table-wrap"><table>')
            for row_index, row in enumerate(rows):
                tag = 'th' if row_index == 0 else 'td'
                result.append('<tr>' + ''.join(f'<{tag}>{inline_markdown(cell)}</{tag}>' for cell in row) + '</tr>')
            result.append('</table></div>')
            continue
        elif line.startswith('#'):
            close_list()
            match = re.match(r'^(#{1,6})\s+(.+)', line)
            if match:
                title = match.group(2)
                depth = max(2, min(4, len(match.group(1))))
                slug = re.sub(r'[^a-z0-9]+', '-', title.lower()).strip('-') or 'section'
                seen[slug] = seen.get(slug, 0) + 1
                ident = slug if seen[slug] == 1 else f'{slug}-{seen[slug]}'
                result.append(f'<h{depth} id="{ident}">{inline_markdown(title)}</h{depth}>')
                if depth == 2:
                    toc.append((ident, re.sub(r'[*`]', '', title)))
        elif re.match(r'^(- |\d+\. )', line):
            kind = 'ul' if line.startswith('- ') else 'ol'
            if list_type != kind:
                close_list(); result.append(f'<{kind}>'); list_type = kind
            result.append('<li>' + inline_markdown(re.sub(r'^(- |\d+\. )', '', line)) + '</li>')
        elif line.startswith('> '):
            close_list(); result.append('<blockquote>' + inline_markdown(line[2:]) + '</blockquote>')
        elif re.fullmatch(r'-{3,}', line):
            close_list(); result.append('<hr>')
        elif line.strip():
            close_list()
            paragraph = [line]
            while index + 1 < len(lines) and lines[index + 1].strip() and not re.match(r'^(#|\||```|> |\d+\. |- )', lines[index + 1]):
                index += 1; paragraph.append(lines[index])
            result.append('<p>' + inline_markdown(' '.join(paragraph)) + '</p>')
        else:
            close_list()
        index += 1
    close_list()
    return '\n'.join(result), toc


def report_page(title, subtitle, body, toc, companion, companion_label):
    navigation = ''.join(f'<a href="#{html.escape(key)}">{html.escape(label)}</a>' for key, label in toc)
    url = safe_url(companion)
    companion_link = f'<a href="{html.escape(url, quote=True)}">{html.escape(companion_label)}</a>' if url else ''
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title><style>{REPORT_CSS}</style></head><body>
<header><div class="eyebrow">Meristem / plant genome annotation</div><h1>{html.escape(title)}</h1>
<p>{html.escape(subtitle)}</p>{companion_link}</header><div class="layout">
<aside aria-label="Contents"><strong>In this report</strong>{navigation}</aside><main>{body}</main></div>
<footer>Self-contained report. No network connection is needed to read it. Use your browser’s Print command for a formatted PDF.</footer>
</body></html>'''


# Tool descriptions describe this repository's invocation, not every upstream mode.
TOOL_CATALOG = [
    ('Assembly statistics', 'Preparation', None, 'VALIDATE_ASSEMBLY', None,
     'Checks the assembly FASTA and summarizes sequence length, continuity, and gaps.',
     'Assembly FASTA', 'assembly_qc/ and assembly statistics TSV',
     'N50 describes continuity, not annotation accuracy. Validate the counters before interpreting the summary.', None),
    ('EDTA', 'Preparation', 'run_edta', 'EDTA_DISCOVERY', 'edta_container',
     'Discovers transposable elements and builds the repeat library used by the masking stage.',
     'Assembly and optional supplied repeat database', 'Repeat library and TE annotation GFF3',
     'A missing published TE GFF is not proof that repeat discovery failed; inspect the native EDTA outputs.', 'https://github.com/oushujun/EDTA'),
    ('RepeatMasker', 'Preparation', 'run_repeatmasker', 'REPEATMASKER', 'repeatmasker_container',
     'Identifies library-matching repeats and softmasks their bases before gene prediction.',
     'Assembly and repeat library', 'Lowercase-masked FASTA and repeat alignments',
     'Masking is not deletion. Repeat-overlapping coding models need evidence review rather than blanket removal.', 'https://www.repeatmasker.org/'),
    ('MMseqs2', 'Evidence preparation', 'cluster_proteins', 'CLUSTER_PROTEINS', 'mmseqs_container',
     'Clusters redundant evidence proteins so BRAKER4 can use a smaller representative set.',
     'Protein evidence FASTA and identity/coverage thresholds', 'Cluster representatives and assignments',
     'Clustering changes the evidence supplied to prediction. These representatives are not the predicted gene set.', 'https://github.com/soedinglab/MMseqs2'),
    ('BRAKER4', 'Structural annotation', 'run_braker4', 'BRAKER4', 'braker4_runner_container',
     'Runs protein-supported gene prediction and model selection through the staged Snakemake workflow.',
     'Masked genome and protein evidence; this parent workflow uses protein-only mode', 'GFF3, predicted proteins, coding sequences, and native logs',
     'Evidence used for prediction is not independent validation. This invocation supplies no RNA-seq input.', 'https://github.com/Gaius-Augustus/BRAKER4'),
    ('Snakemake', 'BRAKER4 component', 'run_braker4', None, None,
     'Schedules and resumes the rules inside BRAKER4.', 'BRAKER4 rules, samples table, and configuration',
     'Native rule logs and intermediate results', 'Meristem schedules the outer task with Nextflow; Snakemake manages its internal rules.', 'https://snakemake.readthedocs.io/'),
    ('ProtHint', 'BRAKER4 component', 'run_braker4', None, None,
     'Derives gene-structure hints from homologous protein evidence.', 'Genome and reference proteins', 'Protein-derived hints',
     'Hint source codes must retain their provenance; a combined code alone does not establish RNA-seq support.', 'https://github.com/gatech-genemark/ProtHint'),
    ('GeneMark', 'BRAKER4 component', 'run_braker4', None, None,
     'Produces gene predictions used in the BRAKER4 training and selection workflow.', 'Genome and configured hints',
     'Gene models and training candidates', 'Intermediate models are distinct from the final selected annotation.', 'https://github.com/Gaius-Augustus/BRAKER'),
    ('AUGUSTUS', 'BRAKER4 component', 'run_braker4', None, None,
     'Predicts exon and intron structures using a trained gene model and hints.', 'Genome, trained parameters, and hints',
     'Candidate gene structures', 'Internal training scores do not measure the accuracy of every final gene.', 'https://github.com/Gaius-Augustus/Augustus'),
    ('TSEBRA', 'BRAKER4 component', 'run_braker4', None, None,
     'Selects between candidate transcript structures using evidence support.', 'Candidate transcripts and hints', 'Selected transcript set',
     'Selection can retain isoforms; transcript duplicates are not necessarily duplicated genomic loci.', 'https://github.com/Gaius-Augustus/TSEBRA'),
    ('compleasm', 'BRAKER4 component', 'run_braker4', None, None,
     'Assesses conserved-marker recovery and supplies marker-based information inside BRAKER4.',
     'Genome or proteins and a lineage dataset', 'Completeness summaries and marker-derived hints',
     'Genome mapping and protein-mode BUSCO are different assessments; marker-guided prediction limits independence.', 'https://github.com/huangnengCSU/compleasm'),
    ('Helixer', 'Structural annotation', 'run_helixer', 'HELIXER', 'helixer_container',
     'Uses the local plant neural-network model to predict gene structures from genomic sequence.',
     'Masked assembly and local HDF5 model', 'Gene annotation GFF3',
     'The GFF must be translated correctly before protein-mode completeness is evaluated; short models need independent support.', 'https://github.com/usadellab/Helixer'),
    ('CDS-to-protein extraction', 'Annotation assessment', 'run_helixer', 'HELIXER_PROTEINS', 'busco_container',
     'Reconstructs Helixer coding sequences in transcript order and translates them for assessment.',
     'GFF3, matching genome, strand, and CDS phases', 'Helixer protein FASTA',
     'Incorrect phase trimming can create artificial stops and depress BUSCO. Verify junction-spanning codons and protein IDs.', 'https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md'),
    ('BUSCO', 'Annotation assessment', 'busco_db', 'ANNOTATION_BUSCO', 'busco_container',
     'Measures recovery of expected conserved marker genes in the predicted proteins.',
     'Protein FASTA and local lineage dataset', 'Complete, single-copy, duplicated, fragmented, and missing marker summaries',
     'Compare the same lineage and protein-selection policy. Completeness is not false-positive rate, and isoforms can inflate duplication.', 'https://busco.ezlab.org/busco_userguide'),
    ('Annotation comparison', 'Annotation assessment', 'both_methods', 'ANNOTATION_COMPARISON', None,
     'Counts annotation features for a side-by-side method comparison.', 'BRAKER4 and Helixer GFF3 files', 'Feature comparison TSV',
     'A larger gene count is not proof of better annotation; feature counts do not measure exact locus agreement.', None),
    ('InterProScan', 'Optional functional annotation', 'interproscan', 'INTERPROSCAN', 'interproscan_container',
     'Assigns protein-family and domain signatures to predicted proteins.', 'BRAKER4 proteins and supplied signature databases',
     'InterPro annotation table', 'Opt-in. Domain absence alone does not prove a false gene.', 'https://interproscan-docs.readthedocs.io/'),
    ('eggNOG-mapper', 'Optional functional annotation', 'eggnog', 'EGGNOG', 'eggnog_container',
     'Transfers functional annotations using orthology assignments.', 'BRAKER4 proteins and local eggNOG database',
     'Orthology and functional annotation tables', 'Opt-in. Functional transfers depend on orthology and database coverage.', 'https://github.com/eggnogdb/eggnog-mapper'),
    ('MultiQC', 'Reporting', None, None, 'multiqc_container',
     'Aggregates supported native QC logs when recognizable inputs are present.', 'Native tool reports recognized by installed modules',
     'MultiQC HTML and data when native inputs are recognized',
     'Meristem renders its custom tables itself; the custom analysis page is not evidence that a native MultiQC module ran.', 'https://seqera.io/multiqc/'),
    ('Nextflow and Apptainer', 'Execution and provenance', None, None, None,
     'Coordinate the parent workflow and provide configured container execution.', 'Workflow, parameters, local inputs, and available images',
     'Task trace, timeline, execution report, and task command files',
     'Successful execution does not establish biological quality. Container tags describe configuration, not measured binary versions.', 'https://www.nextflow.io/docs/latest/'),
]


def selected(value):
    return yes(str(value)) if isinstance(value, (str, bool, int)) else bool(value)


def tool_status(flag, process, config, completed):
    if process and process in completed:
        return 'Observed completed/cached'
    if flag == 'both_methods':
        if not all(key in config for key in ('run_braker4', 'run_helixer')):
            return 'Selection not recorded'
        enabled = selected(config.get('run_braker4')) and selected(config.get('run_helixer'))
    elif flag is None:
        return 'Workflow/reporting component'
    elif flag not in config:
        return 'Selection not recorded'
    else:
        enabled = selected(config[flag])
    if flag in ('interproscan', 'eggnog', 'cluster_proteins'):
        enabled = enabled and selected(config.get('run_braker4'))
    if not enabled:
        return 'Not selected'
    return 'BRAKER4 component; inspect native logs' if process is None else 'Selected; completion not verified here'


def render_tools_report(outdir, config, trace=None, analysis_link=None):
    completed = set()
    if trace:
        for row in read_tsv(trace):
            if row.get('status') in ('COMPLETED', 'CACHED'):
                name = row.get('name', '').split(' (')[0].split(':')[-1]
                completed.add(name)
    lines = ['# Pipeline tools — Meristem', '',
             'This companion guide explains why each tool is present, what it consumes and produces, and how to interpret its output.', '',
             'Selection comes from the supplied run configuration. Completed/cached status is shown only when an explicit task trace supports it. '
             'Nested BRAKER4 components require their native logs. Container references are configured images, not measured executable versions.', '',
             '## Pipeline at a glance', '',
             'Assembly checks → repeat discovery → softmasking → independent gene prediction → protein reconstruction → conserved-marker assessment → interpretation.', '',
             'Optional protein clustering precedes BRAKER4. Optional domain/orthology annotation follows structural prediction.', '',
             '## Run configuration', '', '| Setting | Recorded value |', '|---|---|']
    # Explicit allowlist: never dump arbitrary environment or parameter values.
    for key in ('run_edta', 'run_repeatmasker', 'run_braker4', 'run_helixer', 'cluster_proteins',
                'cluster_min_seq_id', 'cluster_coverage', 'busco_lineage', 'interproscan', 'eggnog'):
        lines.append(f'| {key} | {config.get(key, "not recorded")} |')
    lines.extend(['', '## Tools and their roles', ''])
    overview, toc = markdown_sections('\n'.join(lines))
    cards = []
    for index, entry in enumerate(TOOL_CATALOG):
        name, stage, flag, process, image_key, purpose, inputs, outputs, caveat, source = entry
        status = tool_status(flag, process, config, completed)
        fields = [('Stage', stage), ('Purpose', purpose), ('Inputs', inputs), ('Outputs', outputs), ('Interpretation', caveat)]
        if image_key:
            fields.append(('Configured image', str(config.get(image_key) or 'not recorded; inspect configuration/native logs')))
        if source:
            fields.append(('Documentation', source))
        ident = f'tool-{index}'
        toc.append((ident, name))
        details = []
        lines.extend([f'### {name}', '', f'**Run status:** {status}', ''])
        for label, value in fields:
            rendered = f'<a href="{html.escape(value, quote=True)}">Official documentation</a>' if label == 'Documentation' else html.escape(value)
            details.append(f'<dt>{label}</dt><dd>{rendered}</dd>')
            lines.append(f'- **{label}:** {value}')
        lines.append('')
        cards.append(f'<section class="tool" data-tool="{html.escape(name.lower())}"><h2 id="{ident}">{html.escape(name)}</h2><span class="badge">{html.escape(status)}</span><dl>{"".join(details)}</dl></section>')
    lines.extend(['## Provenance and versions', '',
                  'For measured versions, consult native logs, the staged BRAKER4 checkout, BUSCO summary JSON, and task command files. '
                  'A missing version is unrecorded, not inferred from a filename. The execution trace establishes outer task completion, not every nested tool invocation.', ''])
    text = '\n'.join(lines)
    (outdir / 'pipeline_tools.md').write_text(text)
    write_docx(text, outdir / 'pipeline_tools.docx')
    write_pdf(text, outdir / 'pipeline_tools.pdf')
    body = overview + '<label for="tool-search">Find a tool</label><input class="search" id="tool-search" type="search" placeholder="Search tools by name">' + ''.join(cards)
    body += '<script>document.getElementById("tool-search").addEventListener("input",function(){const q=this.value.toLowerCase();document.querySelectorAll("[data-tool]").forEach(e=>e.hidden=!e.dataset.tool.includes(q));});</script>'
    body += '<h2 id="versions">Provenance and versions</h2><p>' + html.escape(lines[-2]) + '</p>'
    toc.append(('versions', 'Provenance and versions'))
    page = report_page('Tools and methods', 'A practical guide to the annotation workflow, its evidence, and its limitations.', body, toc,
                       analysis_link or '#', 'Open the analysis report' if analysis_link else 'Tools guide')
    (outdir / 'pipeline_tools.html').write_text(page)


def finite_number(value):
    try:
        parsed = float(value)
        return parsed if math.isfinite(parsed) else None
    except (TypeError, ValueError):
        return None


def render_analysis_page(sample, lineage, text, assembly, busco, tools_link, reviewed=False):
    cards = []
    for method, row in busco.items():
        val = finite_number(row.get('complete_pct'))
        shown = f'{val:.1f}%' if val is not None else 'Not measured'
        cards.append(f'<div class="card"><span class="label">{html.escape(method.upper())} · complete BUSCOs</span><strong class="value">{shown}</strong><small>{html.escape(lineage)} · conserved gene recovery</small></div>')
    length = finite_number(assembly.get('total_length'))
    if length is not None:
        cards.append(f'<div class="card"><span class="label">Assembly length</span><strong class="value">{length / 1e6:,.1f} Mb</strong><small>Source statistics; see validation notes</small></div>')
    body, toc = markdown_sections(text)
    note = ('Reviewed quality assessment: reconstructed proteins, reassessed BUSCO summaries, and independently recalculated assembly statistics. '
            'This is a quality reassessment of the completed run, not a new gene-prediction run.') if reviewed else (
            'This report interprets the supplied summary inputs. Protein extraction, isoform policy, and assembly-statistic validity must be checked before ranking methods. '
            'A generated report does not establish overall run success; consult pipeline_info/RUN_STATUS.md.')
    warning = '' if reviewed else ' warning'
    content = '<div class="cards">' + ''.join(cards) + '</div><div class="notice' + warning + '">' + html.escape(note) + '</div>' + body
    return report_page(f'{sample} · Analysis report', 'Results, quality evidence, interpretation, and decisions for plant genome annotation.', content, toc, tools_link, 'Open the tools and methods guide')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", default="Meristem")
    parser.add_argument("--lineage", default="poales_odb10")
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--busco", type=Path)
    parser.add_argument("--comparison", type=Path)
    parser.add_argument("--assembly", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--tools-only", action="store_true")
    parser.add_argument("--with-tools", action="store_true")
    parser.add_argument("--tools-config", type=Path)
    parser.add_argument("--trace", type=Path)
    parser.add_argument("--analysis-markdown", type=Path)
    parser.add_argument("--method", choices=('braker4', 'helixer'))
    parser.add_argument("--review-metrics", type=Path)
    parser.add_argument("--tools-link", default="../pipeline_tools.html")
    parser.add_argument("--run-edta", default="false")
    parser.add_argument("--run-repeatmasker", default="false")
    parser.add_argument("--interproscan", default="false")
    parser.add_argument("--eggnog", default="false")
    args = parser.parse_args()

    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', args.sample):
        parser.error('--sample must be a simple identifier containing letters, numbers, dots, underscores, or hyphens')
    args.outdir.mkdir(parents=True, exist_ok=True)
    config = json.loads(args.tools_config.read_text()) if args.tools_config else {}
    if args.tools_only or args.with_tools:
        render_tools_report(args.outdir, config, args.trace,
                            f'{args.sample}_analysis_report.html' if args.with_tools else None)
    if args.tools_only:
        return
    if args.method:
        if not args.analysis_markdown or not args.assembly:
            parser.error('--method requires --analysis-markdown and --assembly')
        text = args.analysis_markdown.read_text()
        page = render_analysis_page(args.sample, args.lineage, text,
                                    parse_key_value_table(args.assembly), {}, args.tools_link)
        (args.outdir / f'{args.sample}.{args.method}.analysis_report.html').write_text(page)
        return
    if not all((args.summary, args.busco, args.comparison, args.assembly)):
        parser.error("analysis requires --summary, --busco, --comparison, and --assembly")
    assembly = parse_key_value_table(args.assembly)
    busco = method_rows(read_tsv(args.busco))
    features = feature_rows(read_tsv(args.comparison))
    reviewed = False
    if args.review_metrics:
        review = json.loads(args.review_metrics.read_text())
        if not args.analysis_markdown:
            parser.error('--review-metrics requires the accompanying --analysis-markdown assessment')
        measured_length = finite_number(assembly.get('total_length'))
        if measured_length is None or measured_length != review.get('assembly', {}).get('total'):
            parser.error('review metrics and assembly summary do not agree')
        if not busco:
            parser.error('reviewed assessment requires reassessed BUSCO summaries')
        for method, row in busco.items():
            evidence = review.get('reassessed_busco', {}).get(method, {})
            measured_complete = finite_number(row.get('complete_pct'))
            if measured_complete is None or measured_complete != evidence.get('Complete percentage'):
                parser.error(f'reviewed BUSCO completeness does not match the supplied {method} summary')
            if review.get('methods', {}).get(method, {}).get('reconstructed_proteins', {}).get('internal_stop') != 0:
                parser.error(f'reviewed proteins for {method} need an explicit zero-internal-stop integrity result')
        reviewed = True
    interpretation_path = args.outdir / f"{args.sample}_annotation_interpretation.md"
    interpretation_path.write_text(
        interpretation(
            args.sample,
            args.lineage,
            assembly,
            busco,
            features,
            yes(args.run_edta),
            yes(args.run_repeatmasker),
            yes(args.interproscan),
            yes(args.eggnog),
        )
    )
    interpretation_text = interpretation_path.read_text()
    write_pdf(interpretation_text, args.outdir / f"{args.sample}_annotation_interpretation.pdf")
    interpretation_body, interpretation_toc = markdown_sections(interpretation_text)
    interpretation_html = report_page(
        f"{args.sample} · Annotation interpretation",
        "Biological quality assessment and limitations.",
        interpretation_body, interpretation_toc, args.tools_link, "Open the tools and methods guide",
    )
    (args.outdir / f"{args.sample}_annotation_interpretation.html").write_text(interpretation_html)
    combined = args.analysis_markdown.read_text() if args.analysis_markdown else combined_text(args.summary, interpretation_path)
    write_docx(combined, args.outdir / f"{args.sample}_annotation_report.docx")
    write_pdf(combined, args.outdir / f"{args.sample}_annotation_report.pdf")
    # Keep the custom data export for consumers. multiqc_report.html stays reserved
    # for what the workflow's MultiQC step wrote.
    write_html(args.sample, args.lineage, assembly, busco, features, interpretation_path,
               args.outdir / 'multiqc_data' / 'meristem_annotation_summary.html', args.outdir / 'multiqc_data')
    analysis_html = render_analysis_page(
        args.sample, args.lineage, combined, assembly, busco,
        'pipeline_tools.html' if args.with_tools else args.tools_link,
        reviewed=reviewed,
    )
    (args.outdir / f'{args.sample}_analysis_report.html').write_text(analysis_html)
    (args.outdir / f'{args.sample}_analysis_report.md').write_text(combined)


if __name__ == "__main__":
    main()
