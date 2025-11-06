"""
Sanger Sequencing Alignment Tool (Streamlit)

Vibe coded by Tijmen van Slageren 06/11/2025

Uploads a Sanger read (SEQ result from Genewiz) and a reference sequence (.dna from Snapgene).
Performs pairwise alignment, extracts SNPs/INS/DELs and shows a highlighted alignment
and a summary table.
"""

import io
import html
import streamlit as st
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import streamlit.components.v1 as components


st.title("Site Directed Mutagenesis Sequence Alignment Tool")
st.write("Upload your Sanger sequencing read (seq) and a reference snapgene.dna to identify mutations.")
st.caption("Note: reads are automatically trimmed to first 1000 bases. Option to trim first 30 bases and last 150 bases before alignment. Minimum base quality is fixed at Phred=40 internally. Alignment achieved using Biopython's pairwise2.")

# Put file uploaders in the left-hand sidebar for easier access
# Put file uploaders back in the sidebar (left-hand column)
st.sidebar.header("Inputs")
uploaded_sanger = st.sidebar.file_uploader("Upload Sanger sequencing read", type=["seq"])
uploaded_ref = st.sidebar.file_uploader("Upload reference sequence (SnapGene .dna)", type=["dna", "snapgene"])

# Minimum base quality (Phred) used internally for trimming and variant filtering.
# Fixed at 40 as requested (no UI control).
qc_cutoff = 40
trim_enable = st.checkbox("Trim low-quality ends (fixed 30/150)", value=True)
flank_size = st.number_input("Flank bases to show around each mutation", min_value=0, max_value=200, value=42)


def read_uploaded_record(uploaded_file):
    """Read an uploaded Streamlit file (UploadedFile) into a SeqRecord.
    Tries common chromatogram (abi, scf) and text formats (fasta, fastq...).
    Returns SeqRecord or raises ValueError if nothing readable found.
    """
    data = uploaded_file.getvalue()
    # Try binary chromatogram formats
    bio = io.BytesIO(data)
    for fmt in ("abi", "ab1", "scf"):
        try:
            bio.seek(0)
            rec = SeqIO.read(bio, fmt)
            if rec and len(rec.seq) > 0:
                return rec
        except Exception:
            pass

    # Try text formats (decode)
    try:
        text = data.decode("utf-8", errors="replace")
    except Exception:
        text = None

    if text is not None:
        for fmt in ("fasta", "fastq", "genbank", "embl"):
            try:
                sio = io.StringIO(text)
                rec = SeqIO.read(sio, fmt)
                if rec and len(rec.seq) > 0:
                    return rec
            except Exception:
                pass

        # fallback: plain raw sequence text (no header)
        seq_text = "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))
        if seq_text:
            return SeqRecord(Seq(seq_text), id=uploaded_file.name)

    raise ValueError("Could not parse uploaded sequence file")


def align_best(ref_seq, query_seq):
    """Return best global alignment (ref_aligned, query_aligned, score) using PairwiseAligner.

    Uses the same scoring as the previous pairwise2.globalms call: match=2, mismatch=-1,
    gap open=-5, gap extend=-0.5. Returns aligned reference and query strings with '-' for gaps
    plus the alignment score. If no alignment is found, returns None.
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5

    alns = aligner.align(ref_seq, query_seq)
    if len(alns) == 0:
        return None
    aln = alns[0]
    score = aln.score

    # Reconstruct gapped sequences from aligned blocks
    aligned_blocks = aln.aligned
    ref_blocks = aligned_blocks[0]
    qry_blocks = aligned_blocks[1]

    ref_idx = 0
    qry_idx = 0
    aln_ref = []
    aln_qry = []

    for (r_start, r_end), (q_start, q_end) in zip(ref_blocks, qry_blocks):
        # add any gaps/unaligned region before the block
        while ref_idx < r_start or qry_idx < q_start:
            if ref_idx < r_start and qry_idx < q_start:
                # both advanced - take bases (this shouldn't usually happen between blocks)
                aln_ref.append(ref_seq[ref_idx])
                aln_qry.append(query_seq[qry_idx])
                ref_idx += 1
                qry_idx += 1
            elif ref_idx < r_start:
                aln_ref.append(ref_seq[ref_idx])
                aln_qry.append('-')
                ref_idx += 1
            else:
                aln_ref.append('-')
                aln_qry.append(query_seq[qry_idx])
                qry_idx += 1

        # add the aligned block (no gaps inside)
        for i in range(r_end - r_start):
            aln_ref.append(ref_seq[r_start + i])
            aln_qry.append(query_seq[q_start + i])
        ref_idx = r_end
        qry_idx = q_end

    # append any trailing tails
    while ref_idx < len(ref_seq) or qry_idx < len(query_seq):
        if ref_idx < len(ref_seq) and qry_idx < len(query_seq):
            aln_ref.append(ref_seq[ref_idx])
            aln_qry.append(query_seq[qry_idx])
            ref_idx += 1
            qry_idx += 1
        elif ref_idx < len(ref_seq):
            aln_ref.append(ref_seq[ref_idx])
            aln_qry.append('-')
            ref_idx += 1
        else:
            aln_ref.append('-')
            aln_qry.append(query_seq[qry_idx])
            qry_idx += 1

    return ''.join(aln_ref), ''.join(aln_qry), score


def parse_variants(aln_ref, aln_query, phred_quals=None):
    variants = []
    i = 0
    ref_pos = 0
    query_pos = 0

    while i < len(aln_ref):
        r = aln_ref[i]
        q = aln_query[i]

        if r != "-" and q != "-":
            ref_pos += 1
            query_pos += 1
            if r != q:
                qual = None
                if phred_quals and (query_pos - 1) < len(phred_quals):
                    qual = phred_quals[query_pos - 1]
                variants.append({
                    "type": "SNP",
                    "ref_pos": ref_pos,
                    "ref": r,
                    "query": q,
                    "query_pos": query_pos,
                    "qual": qual,
                })
            i += 1
            continue

        if r == "-" and q != "-":
            ins_bases = []
            ins_quals = []
            while i < len(aln_ref) and aln_ref[i] == "-" and aln_query[i] != "-":
                ins_bases.append(aln_query[i])
                if phred_quals and query_pos < len(phred_quals):
                    ins_quals.append(phred_quals[query_pos])
                query_pos += 1
                i += 1
            variants.append({
                "type": "INS",
                "ref_pos": ref_pos,
                "ref": "-",
                "query": "".join(ins_bases),
                "query_pos": query_pos - len(ins_bases) + 1,
                "qual": int(sum(ins_quals)/len(ins_quals)) if ins_quals else None,
            })
            continue

        if r != "-" and q == "-":
            del_bases = []
            while i < len(aln_ref) and aln_ref[i] != "-" and aln_query[i] == "-":
                del_bases.append(aln_ref[i])
                ref_pos += 1
                i += 1
            variants.append({
                "type": "DEL",
                "ref_pos": ref_pos - len(del_bases) + 1,
                "ref": "".join(del_bases),
                "query": "-",
                "query_pos": None,
                "qual": None,
            })
            continue

    return variants


def highlight_html(aln_ref, aln_query):
    out_ref = []
    out_query = []
    for r, q in zip(aln_ref, aln_query):
        if r == "-" and q != "-":
            out_ref.append('<span style="color:#bbb">-</span>')
            out_query.append(f'<span style="background:#cfefff;color:#000">{html.escape(q)}</span>')
        elif r != "-" and q == "-":
            out_ref.append(f'<span style="background:#ffd7b5;color:#000">{html.escape(r)}</span>')
            out_query.append('<span style="color:#bbb">-</span>')
        elif r == q:
            out_ref.append(html.escape(r))
            out_query.append(html.escape(q))
        else:
            out_ref.append(f'<span style="background:#ffd6d6;color:#000">{html.escape(r)}</span>')
            out_query.append(f'<span style="background:#ffd6d6;color:#000">{html.escape(q)}</span>')
    html_text = "<pre style='font-family:monospace'>REF:  " + "".join(out_ref) + "\nQUERY: " + "".join(out_query) + "</pre>"
    return html_text


def format_alignment_html(aln_ref, aln_query, width=80, ref_pos_map=None, query_pos_map=None):
    """Produce a block-formatted HTML alignment with coordinates and coloring.
    Shows reference and query in fixed-width font, with a middle match line and a legend.
    """
    # compute position maps for ref and query (1-based positions or '-') if not provided
    if ref_pos_map is None or query_pos_map is None:
        ref_pos_map = []
        query_pos_map = []
        ref_ctr = 0
        query_ctr = 0
        for r, q in zip(aln_ref, aln_query):
            if r != "-":
                ref_ctr += 1
                ref_pos_map.append(str(ref_ctr))
            else:
                ref_pos_map.append("-")
            if q != "-":
                query_ctr += 1
                query_pos_map.append(str(query_ctr))
            else:
                query_pos_map.append("-")

    # helper for coloring bases
    def color_base(r, q):
        if r == "-" and q != "-":
            return (f'<span style="background:#cfefff;color:#000">{html.escape(q)}</span>', '<span style="color:#bbb">-</span>')
        if r != "-" and q == "-":
            return (f'<span style="background:#ffd7b5;color:#000">{html.escape(r)}</span>', '<span style="color:#bbb">-</span>')
        if r == q:
            return (html.escape(r), html.escape(q))
        return (f'<span style="background:#ffd6d6;color:#000">{html.escape(r)}</span>', f'<span style="background:#ffd6d6;color:#000">{html.escape(q)}</span>')

    blocks = []
    length = len(aln_ref)
    for start in range(0, length, width):
        end = min(start + width, length)
        ref_bases_html = []
        query_bases_html = []
        match_line = []
        ref_start_pos = None
        ref_end_pos = None
        query_start_pos = None
        query_end_pos = None

        for i in range(start, end):
            r = aln_ref[i]
            q = aln_query[i]
            r_html, q_html = color_base(r, q)
            ref_bases_html.append(r_html)
            query_bases_html.append(q_html)
            if r == q and r != "-":
                match_line.append("|")
            else:
                match_line.append(" ")
            # positions
            rp = ref_pos_map[i]
            qp = query_pos_map[i]
            if rp != "-":
                if ref_start_pos is None:
                    ref_start_pos = int(rp)
                ref_end_pos = int(rp)
            if qp != "-":
                if query_start_pos is None:
                    query_start_pos = int(qp)
                query_end_pos = int(qp)

        # format position labels
        ref_label = f"{ref_start_pos if ref_start_pos is not None else '-'}-{ref_end_pos if ref_end_pos is not None else '-'}"
        query_label = f"{query_start_pos if query_start_pos is not None else '-'}-{query_end_pos if query_end_pos is not None else '-'}"

        block_html = (
            f"<div style='font-family:monospace; margin-bottom:8px;'>"
            f"<div><strong>REF {ref_label:>12}</strong></div>"
            f"<div style='white-space:pre;'>{''.join(ref_bases_html)}</div>"
            f"<div style='white-space:pre;color:#666'>{''.join(match_line)}</div>"
            f"<div><strong>QRY {query_label:>12}</strong></div>"
            f"<div style='white-space:pre;'>{''.join(query_bases_html)}</div>"
            f"</div>"
        )
        blocks.append(block_html)

    legend = ("<div style='font-family:monospace; margin-top:6px;'>"
              "<strong>Legend:</strong> <span style='background:#ffd6d6'>SNP</span> "
              "<span style='background:#cfefff'>INS</span> <span style='background:#ffd7b5'>DEL</span>"
              "</div>")

    full_html = """
    <div>
    <style> .align-block { max-width: 100%; } </style>
    <div class='align-block'>
    """ + "".join(blocks) + "</div></div>"
    # include legend by default at the end (can be suppressed by caller by not using it)
    return full_html, legend
    


def compute_pos_maps(aln_ref, aln_query):
    """Return ref_pos_map and query_pos_map lists (strings or '-') for an alignment."""
    ref_pos_map = []
    query_pos_map = []
    ref_ctr = 0
    query_ctr = 0
    for r, q in zip(aln_ref, aln_query):
        if r != "-":
            ref_ctr += 1
            ref_pos_map.append(str(ref_ctr))
        else:
            ref_pos_map.append("-")
        if q != "-":
            query_ctr += 1
            query_pos_map.append(str(query_ctr))
        else:
            query_pos_map.append("-")
    return ref_pos_map, query_pos_map


def find_alignment_indices_for_variants(variants, ref_pos_map, query_pos_map):
    """Given filtered variants (with ref_pos and/or query_pos), return set of alignment indices containing them."""
    indices = set()
    for v in variants:
        refp = v.get("ref_pos")
        queryp = v.get("query_pos")
        for i, (rp, qp) in enumerate(zip(ref_pos_map, query_pos_map)):
            if rp != "-" and refp is not None and int(rp) == int(refp):
                indices.add(i)
            if qp != "-" and queryp is not None and int(qp) == int(queryp):
                indices.add(i)
    return sorted(indices)


def build_mutation_region_html(aln_ref, aln_query, variant_indices, ref_pos_map, query_pos_map, flank=10):
    """Build HTML showing only regions around the given variant alignment indices.
    Merges overlapping windows and returns combined HTML (using format_alignment_html on slices).
    """
    if not variant_indices:
        return ""
    length = len(aln_ref)
    # build intervals
    intervals = []
    for idx in variant_indices:
        a = max(0, idx - flank)
        b = min(length, idx + flank + 1)
        intervals.append((a, b))
    # merge
    intervals.sort()
    merged = []
    cur_a, cur_b = intervals[0]
    for a, b in intervals[1:]:
        if a <= cur_b:
            cur_b = max(cur_b, b)
        else:
            merged.append((cur_a, cur_b))
            cur_a, cur_b = a, b
    merged.append((cur_a, cur_b))

    parts = []
    legend = ""
    for a, b in merged:
        sub_ref = aln_ref[a:b]
        sub_qry = aln_query[a:b]
        sub_ref_map = ref_pos_map[a:b]
        sub_qry_map = query_pos_map[a:b]
        html_part, legend_part = format_alignment_html(sub_ref, sub_qry, width=b - a, ref_pos_map=sub_ref_map, query_pos_map=sub_qry_map)
        parts.append(html_part)
        # capture legend from last call (they're identical)
        legend = legend_part
    # return blocks and legend separately so caller can render legend once and make regions scrollable
    blocks_html = "<div>" + "".join(parts) + "</div>"
    return blocks_html, legend


def compute_max_windows(aln_ref, aln_query, ref_pos_map, query_pos_map, max_mismatches=5):
    """Compute windows (start,end) in alignment with at most max_mismatches mismatches.
    Returns a list of dicts with alignment indices, lengths and corresponding ref/query coords.
    Uses a two-pointer approach to find the maximal right for each left.
    """
    N = len(aln_ref)
    results = []
    r = 0
    mismatch = 0

    def is_mismatch(i):
        a = aln_ref[i]
        b = aln_query[i]
        return (a != "-" and b != "-" and a != b)

    for l in range(N):
        # advance r as far as possible while keeping mismatches <= max_mismatches
        while r < N and (mismatch + (1 if is_mismatch(r) else 0)) <= max_mismatches:
            if is_mismatch(r):
                mismatch += 1
            r += 1

        # window [l, r-1] is maximal for this left
        if r - 1 >= l:
            # compute mismatch count in window
            mismatches_in_window = 0
            for i in range(l, r):
                if is_mismatch(i):
                    mismatches_in_window += 1

            # find ref/query start/end coords
            ref_start = None
            ref_end = None
            qry_start = None
            qry_end = None
            for i in range(l, r):
                rp = ref_pos_map[i]
                qp = query_pos_map[i]
                if rp != "-":
                    if ref_start is None:
                        ref_start = int(rp)
                    ref_end = int(rp)
                if qp != "-":
                    if qry_start is None:
                        qry_start = int(qp)
                    qry_end = int(qp)

            results.append({
                "aln_start": l,
                "aln_end": r - 1,
                "length": r - l,
                "mismatches": mismatches_in_window,
                "ref_start": ref_start,
                "ref_end": ref_end,
                "qry_start": qry_start,
                "qry_end": qry_end,
            })

        # move left pointer forward: if left was a mismatch, decrement mismatch count
        if r <= l:
            r = l + 1
        else:
            if is_mismatch(l):
                mismatch -= 1

    # filter out degenerate entries with length<=0
    results = [r for r in results if r["length"] > 0]
    # collapse duplicates by keeping the longest window for a given aln_start
    best_by_start = {}
    for item in results:
        s = item["aln_start"]
        if s not in best_by_start or item["length"] > best_by_start[s]["length"]:
            best_by_start[s] = item
    results = list(best_by_start.values())
    # sort by length desc
    results.sort(key=lambda x: x["length"], reverse=True)
    return results


if uploaded_sanger and uploaded_ref:
    # Try reading the uploaded reference as SnapGene (.dna) first, then fallback to FASTA
    ref_rec = None
    ref_bytes = uploaded_ref.getvalue()
    try:
        bio = io.BytesIO(ref_bytes)
        bio.seek(0)
        ref_rec = SeqIO.read(bio, "snapgene")
    except Exception:
        try:
            text = ref_bytes.decode("utf-8", errors="replace")
            ref_rec = SeqIO.read(io.StringIO(text), "fasta")
        except Exception:
            st.error("Could not read reference file as SnapGene .dna or FASTA. Ensure the file is valid.")
            ref_rec = None

    if ref_rec is None:
        st.stop()

    try:
        read_rec = read_uploaded_record(uploaded_sanger)
    except Exception as e:
        st.error(f"Could not parse uploaded Sanger read: {e}")
        st.stop()

    # Hard-limit read to first 1000 bases BEFORE trimming
    max_read_len = 1000
    orig_len = len(read_rec.seq)
    if orig_len > max_read_len:
        read_rec.seq = read_rec.seq[:max_read_len]
        if hasattr(read_rec, "letter_annotations") and "phred_quality" in read_rec.letter_annotations:
            read_rec.letter_annotations["phred_quality"] = read_rec.letter_annotations["phred_quality"][:max_read_len]
        st.write(f"Truncated read to first {max_read_len} bases (original length {orig_len}) before trimming.")

    # perform optional trimming based on phred qualities if available
    def trim_low_quality_ends(rec):
        """Trim fixed low-quality ends from a SeqRecord: first 30 bases and last 150 bases.
        Returns a trimmed SeqRecord and tuple (left_trim, right_trim) counts.
        """
        seq = str(rec.seq)
        n = len(seq)
        if n == 0:
            return rec, (0, 0)

        # fixed trimming amounts
        left_fixed = 30
        right_fixed = 150

        left_base = min(left_fixed, n)
        right_base = max(0, n - right_fixed)

        # ensure we don't trim everything; keep at least one base
        if left_base >= right_base:
            left_base = 0
            right_base = min(n, 1)

        new_seq = seq[left_base:right_base]
        new_rec = SeqRecord(Seq(new_seq), id=rec.id)

        # copy phred qualities for trimmed region if present
        phreds = rec.letter_annotations.get("phred_quality") if hasattr(rec, "letter_annotations") else None
        if phreds and len(phreds) == n:
            new_phreds = phreds[left_base:right_base]
            if new_phreds:
                new_rec.letter_annotations["phred_quality"] = new_phreds

        return new_rec, (left_base, n - right_base)

    if trim_enable:
        read_rec_trimmed, (left_trim, right_trim) = trim_low_quality_ends(read_rec)
        if left_trim or right_trim:
            st.write(f"Trimmed {left_trim} bases from left and {right_trim} bases from right due to low quality.")
        read_rec = read_rec_trimmed

    # reverse-complement functionality removed — reads are used as provided

    

    phreds = read_rec.letter_annotations.get("phred_quality") if hasattr(read_rec, "letter_annotations") else None

    # Align (ref first so aln_ref corresponds to reference)
    # Normalize case to avoid mismatches caused by lowercase letters in reference
    # (SnapGene sometimes stores features in lowercase). Uppercase does not change
    # sequence length, so it's safe for alignment and position mapping.
    aln = align_best(str(ref_rec.seq).upper(), str(read_rec.seq).upper())
    if not aln:
        st.error("No alignment found between reference and read.")
        st.stop()

    aln_ref, aln_query, score = aln

    # plain-text alignment for download (full alignment)
    def build_plain_alignment(aln_ref, aln_query, width=80):
        lines = []
        length = len(aln_ref)
        for start in range(0, length, width):
            end = min(start + width, length)
            ref_part = ''.join(aln_ref[start:end])
            qry_part = ''.join(aln_query[start:end])
            match = ''.join('|' if (r == q and r != '-') else ' ' for r, q in zip(ref_part, qry_part))
            lines.append(f"REF  {start+1}-{end}: {ref_part}")
            lines.append(f"      {match}")
            lines.append(f"QRY  {start+1}-{end}: {qry_part}")
            lines.append("")
        return "\n".join(lines)

    plain_align = build_plain_alignment(aln_ref, aln_query, width=80)
    st.download_button("Download full plain alignment", data=plain_align.encode('utf-8'), file_name='alignment.txt', mime='text/plain')

    variants = parse_variants(aln_ref, aln_query, phred_quals=phreds)

    # Filter by quality threshold if provided (only for SNPs and INS where qual exists)
    filtered = []
    for v in variants:
        if v.get("qual") is None:
            filtered.append(v)
        else:
            if v["qual"] >= qc_cutoff:
                filtered.append(v)

    if filtered:
        # NOTE: 'Identified variants' table and CSV download removed by user request.
        # We still compute mutated regions and show alignment slices around them.
        ref_pos_map, query_pos_map = compute_pos_maps(aln_ref, aln_query)
        variant_indices = find_alignment_indices_for_variants(filtered, ref_pos_map, query_pos_map)
        if not variant_indices:
            st.info("No alignment indices found for the filtered variants — showing full alignment instead.")
            full_html, legend = format_alignment_html(aln_ref, aln_query, width=80, ref_pos_map=ref_pos_map, query_pos_map=query_pos_map)
            # show legend once, then full alignment in a scrollable container
            st.markdown(legend, unsafe_allow_html=True)
            scroll_html = f"<div style='max-height:400px; overflow:auto; border:1px solid #ddd; padding:6px;'>{full_html}</div>"
            components.html(scroll_html, height=420)
        else:
            blocks_html, legend = build_mutation_region_html(aln_ref, aln_query, variant_indices, ref_pos_map, query_pos_map, flank=int(flank_size))
            st.subheader("Alignment regions around mutations")
            # show legend once above scrollable regions
            st.markdown(legend, unsafe_allow_html=True)
            scroll_html = f"<div style='max-height:400px; overflow:auto; border:1px solid #ddd; padding:6px;'>{blocks_html}</div>"
            components.html(scroll_html, height=420)
            # Build a minimal table listing SNP locations (ref and query positions)
            # and the current Sanger read length (post truncation/trimming). This
            # shows only SNPs found in the alignment.
            snp_variants = [v for v in filtered if v.get("type") == "SNP"]
            if snp_variants:
                sanger_length = len(read_rec.seq)
                tbl_rows = []
                for v in snp_variants:
                    ref_base = v.get("ref") if v.get("ref") is not None else "-"
                    qry_base = v.get("query") if v.get("query") is not None else "-"
                    tbl_rows.append({
                        "ref_pos": v.get("ref_pos"),
                        "qry_pos": v.get("query_pos"),
                        "change": f"{ref_base}>{qry_base}",
                        "sanger_length": sanger_length,
                    })
                st.subheader("SNP locations and Sanger read length")
                st.table(pd.DataFrame(tbl_rows))
            else:
                st.write("No SNPs were found in the filtered variants.")
    else:
        st.write("No variants passed the filters.")


# Acknowledgements / citations
st.markdown(
    "---\n"
    "**Acknowledgement:** This app uses Biopython for sequence parsing and alignment. "
    "Please cite: Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., "
    "Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. L. (2009). Biopython: freely available Python tools "
    "for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423. "
    "doi:10.1093/bioinformatics/btp163 — https://biopython.org/"
)

