import pysam
import argparse
from collections import Counter
import sw

# Parameters for clustering and filtering
MIN_THRESHOLD = 50
MAX_THRESHOLD = 200
MIN_NUM_CHARS = 10  # Reduced for more consensus sequences
MIN_CHAR_DENSITY = 0.6  # Reduced for more consensus sequences
MIN_SEQ_LEN = 10  # Minimum length for consensus sequences
MAX_CIRCLE_SIZE = 50000  # Maximum microDNA circle size (bp)
MIN_CIRCLE_SIZE = 50  # Minimum microDNA circle size (bp)
MIN_MAPQ = 20  # Minimum mapping quality for reads

def most_common_char(char_list):
    if not char_list:
        return None
    count = Counter(char_list)
    return count.most_common(1)[0][0]

def count_matching_chars(char, char_list):
    return sum(1 for c in char_list if c == char)

def get_consensus(seq_list):
    consensus = {}
    for seq in seq_list:
        for i in range(len(seq)):
            if i not in consensus:
                consensus[i] = []
            consensus[i].append(seq[i])
    S = ''
    for i in sorted(consensus.keys()):
        mcc = most_common_char(consensus[i])
        mcc_count = count_matching_chars(mcc, consensus[i])
        density = mcc_count / len(consensus[i])
        if mcc_count >= MIN_NUM_CHARS and density >= MIN_CHAR_DENSITY:
            S += mcc
    return S

def get_clipped_seq(read):
    clipped_seq = ''
    pos = 0
    for op, length in read.cigartuples:
        if op == 4:  # Soft clip
            clipped_seq += read.query_sequence[pos:pos+length]
            pos += length
        elif op == 0:  # Match
            pos += length
    return clipped_seq

def identify_breakpoints(bam_file, region_start=None, region_end=None):
    start_clusters = []
    end_clusters = []
    current_start_window = []
    current_end_window = []

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Use fetch to limit to a specific region if provided
        iterator = bam.fetch("NC_000001.10", region_start, region_end) if region_start and region_end else bam
        for read in iterator:
            if read.is_unmapped or read.mapping_quality < MIN_MAPQ:
                continue
            if read.cigartuples and read.cigartuples[0][0] == 4:  # Soft clip at start
                if not current_start_window:
                    current_start_window.append((read.pos, get_clipped_seq(read)))
                elif read.pos == current_start_window[0][0]:
                    current_start_window.append((read.pos, get_clipped_seq(read)))
                else:
                    if MIN_THRESHOLD < len(current_start_window) < MAX_THRESHOLD:
                        pos = current_start_window[0][0]
                        seqs = [s[1][::-1] for s in current_start_window]
                        consensus = get_consensus(seqs)[::-1]
                        if consensus and len(consensus) >= MIN_SEQ_LEN:
                            start_clusters.append((pos, consensus))
                    current_start_window = [(read.pos, get_clipped_seq(read))]
            elif read.cigartuples and read.cigartuples[-1][0] == 4:  # Soft clip at end
                circle_end = read.pos + sum(length for op, length in read.cigartuples if op == 0)
                if not current_end_window:
                    current_end_window.append((circle_end, get_clipped_seq(read)))
                elif circle_end == current_end_window[0][0]:
                    current_end_window.append((circle_end, get_clipped_seq(read)))
                else:
                    if MIN_THRESHOLD < len(current_end_window) < MAX_THRESHOLD:
                        pos = current_end_window[0][0]
                        seqs = [s[1] for s in current_end_window]
                        consensus = get_consensus(seqs)
                        if consensus and len(consensus) >= MIN_SEQ_LEN:
                            end_clusters.append((pos, consensus))
                    current_end_window = [(circle_end, get_clipped_seq(read))]

        # Process final windows
        if MIN_THRESHOLD < len(current_start_window) < MAX_THRESHOLD:
            pos = current_start_window[0][0]
            seqs = [s[1][::-1] for s in current_start_window]
            consensus = get_consensus(seqs)[::-1]
            if consensus and len(consensus) >= MIN_SEQ_LEN:
                start_clusters.append((pos, consensus))
        if MIN_THRESHOLD < len(current_end_window) < MAX_THRESHOLD:
            pos = current_end_window[0][0]
            seqs = [s[1] for s in current_end_window]
            consensus = get_consensus(seqs)
            if consensus and len(consensus) >= MIN_SEQ_LEN:
                end_clusters.append((pos, consensus))

    return start_clusters, end_clusters

def align_sequences(start_pos, end_pos, start_seq, end_seq, ref_seq):
    print(f"Aligning: start_pos={start_pos}, end_pos={end_pos}, start_seq={start_seq}, end_seq={end_seq}")

    # Validate inputs
    if not start_seq or not end_seq:
        print("Error: Empty start_seq or end_seq")
        return None
    if start_pos >= end_pos:
        print("Error: Invalid positions")
        return None

    # Slice reference sequence from cached ref_seq
    try:
        seq = ref_seq[start_pos:end_pos]
    except IndexError:
        print("Error: Invalid reference sequence range")
        return None
    if not seq:
        print("Error: No sequence fetched")
        return None

    # Extract subsequences for alignment
    seq_start = seq[:20]  # Increased length for better alignment
    seq_end = seq[-20:]
    print(f"seq_start={seq_start}, seq_end={seq_end}")

    # Align start_seq and end_seq to seq_start and seq_end
    H = sw.sw_fill_matrix(seq_start, end_seq, -1, -1, 2)
    align_A, align_B, score_start = sw.sw_traceback(H, seq_start, end_seq, -1, -1, 2)
    print(f"Start alignment: A={align_A}, B={align_B}, score={score_start}")

    micro_h = align_A
    c_begin = micro_h + end_seq
    c_end = start_seq + micro_h
    print(f"c_begin={c_begin}, c_end={c_end}")

    # Align c_begin to the start of the sequence
    start_seq_ref = seq[:len(c_begin)]
    if not start_seq_ref or not c_begin:
        print("Error: Empty start_seq_ref or c_begin")
        return None
    H = sw.sw_fill_matrix(start_seq_ref, c_begin, -1, -1, 2)
    align_A_begin, align_B_begin, score_begin = sw.sw_traceback(H, start_seq_ref, c_begin, -1, -1, 2)
    print(f"Begin alignment: A={align_A_begin}, B={align_B_begin}, score={score_begin}")

    # Align c_end to the end of the sequence
    end_seq_ref = seq[-len(c_end):]
    if not end_seq_ref or not c_end:
        print("Error: Empty end_seq_ref or c_end")
        return None
    H = sw.sw_fill_matrix(end_seq_ref, c_end, -1, -1, 2)
    align_A_end, align_B_end, score_end = sw.sw_traceback(H, end_seq_ref, c_end, -1, -1, 2)
    print(f"End alignment: A={align_A_end}, B={align_B_end}, score={score_end}")

    return {
        "chrom": "NC_000001.10",
        "start": start_pos,
        "end": end_pos,
        "start_consensus": start_seq,
        "end_consensus": end_seq,
        "ref_sequence": seq,
        "start_alignment": (align_A_begin, align_B_begin, score_begin),
        "end_alignment": (align_A_end, align_B_end, score_end)
    }

def main():
    parser = argparse.ArgumentParser(description="Detect microDNA circles from BAM file")
    parser.add_argument("--bam", type=str, required=True, help="Path to sorted BAM file")
    parser.add_argument("--ref", type=str, required=True, help="Path to reference genome FASTA")
    parser.add_argument("--output", type=str, required=True, help="Output file for results")
    parser.add_argument("--region-start", type=int, help="Start of region to analyze")
    parser.add_argument("--region-end", type=int, help="End of region to analyze")
    args = parser.parse_args()

    # Load and cache reference sequence
    ref_genome = pysam.FastaFile(args.ref)
    ref_seq = ref_genome.fetch("NC_000001.10")
    ref_genome.close()

    # Step 1: Identify breakpoints
    start_clusters, end_clusters = identify_breakpoints(args.bam, args.region_start, args.region_end)
    print(f"Found {len(start_clusters)} start clusters and {len(end_clusters)} end clusters")

    # Step 2: Sort and filter breakpoint pairs
    start_clusters.sort(key=lambda x: x[0])  # Sort by position
    end_clusters.sort(key=lambda x: x[0])
    results = []
    for i, (start_pos, start_seq) in enumerate(start_clusters):
        # Only check end clusters within MAX_CIRCLE_SIZE
        for end_pos, end_seq in end_clusters:
            if end_pos <= start_pos:
                continue
            if end_pos - start_pos > MAX_CIRCLE_SIZE:
                break  # Since end_clusters is sorted, no further pairs will be valid
            if MIN_CIRCLE_SIZE <= end_pos - start_pos and len(start_seq) >= MIN_SEQ_LEN and len(end_seq) >= MIN_SEQ_LEN:
                result = align_sequences(start_pos, end_pos, start_seq, end_seq, ref_seq)
                if result and result["start_alignment"][2] > 5 and result["end_alignment"][2] > 5:
                    results.append(result)

    print(f"Processed {len(results)} valid alignments")

    # Step 3: Write results to output file
    with open(args.output, "w") as f:
        f.write("Chrom\tStart\tEnd\tStart_Consensus\tEnd_Consensus\tRef_Sequence\tStart_Alignment\tEnd_Alignment\tStart_Score\tEnd_Score\n")
        for res in results:
            f.write(f"{res['chrom']}\t{res['start']}\t{res['end']}\t"
                    f"{res['start_consensus']}\t{res['end_consensus']}\t{res['ref_sequence']}\t"
                    f"{res['start_alignment'][0]}|{res['start_alignment'][1]}\t"
                    f"{res['end_alignment'][0]}|{res['end_alignment'][1]}\t"
                    f"{res['start_alignment'][2]}\t{res['end_alignment'][2]}\n")

    print(f"Results written to {args.output}")

if __name__ == "__main__":
    main()