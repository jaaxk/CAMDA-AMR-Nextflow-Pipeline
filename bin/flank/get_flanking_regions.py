from Bio import SeqIO
from Bio.Seq import Seq
import sys

def parse_blast_output(blast_file):
    hits = []
    with open(blast_file) as f:
        for line in f:
            cols = line.strip().split('\t')
            qseqid = cols[0]
            sseqid = cols[1]
            sstart = int(cols[6])
            send = int(cols[7])
            strand = cols[8] if len(cols) > 8 else '+'
            hits.append((qseqid, sseqid, sstart, send, strand))
    return hits

def main():
    assembly_file = sys.argv[1]
    seq_length = int(sys.argv[2])
    flank_length = seq_length // 2
    blast_file = sys.argv[3]
    output_file = sys.argv[4]

    seqs = []
    assembly_dict = SeqIO.to_dict(SeqIO.parse(assembly_file, 'fasta'))
    hits = parse_blast_output(blast_file)
    for qseqid, sseqid, sstart, send, strand in hits:
        if sseqid not in assembly_dict:
            continue
        
        contig_seq = assembly_dict[sseqid].seq
        seq_len = len(contig_seq)

        mid = (sstart + send) // 2
        start = max(0, mid - flank_length)
        end = min(seq_len, mid + flank_length)
        region_seq = contig_seq[start:end]
        if len(region_seq) < seq_length/2: #Can't be smaller than half the sequence length
            continue
        seqs.append(region_seq)

    if len(seqs) == 0:
        print("ERROR: No sequences found")
        sys.exit(1)

    with open(output_file, 'w') as f:
        for i, seq in enumerate(seqs):
            f.write('>seq_{}\n{}\n'.format(i, seq))

if __name__ == '__main__':
    main()

    
        
        