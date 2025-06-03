import pandas as pd
import gzip

"""
Helper Functions
"""
###convert nucleotide sequence to amino acid sequence
def translate_sequence(sequence):
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    amino_acids = []
    codons = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]

    for codon in codons:
        amino_acid = genetic_code.get(codon.upper(), '')
        amino_acids.append(amino_acid)

    return ''.join(amino_acids)

####Find CDS start and stop in APPRIS FASTA file
def find_cds(transcript, fastaf):
    file_start = fastaf.find(transcript)
    if file_start == -1:
        return []
    file_name = fastaf[file_start:fastaf.find("\\n", file_start)]
    cds_coor = file_name.find("CDS")
    trans = file_name[cds_coor:file_name.find("|", cds_coor)]
    cds = trans[4:].split("-")
    output = [int(cds[0]), int(cds[1])]
    return output

####add value-key pairs to a dictionary
def add_dic(key, value, dict):
    if key in dict.keys():
        dict[key].append(value)
    else:
        dict[key] = [value]






"""
Algorithm to add additional analysis information to ORFs file - nucleotide sequence, 
    amino acid sequence, CDS amino acid sequence, CDS start and stop coordinates
"""
def write_output(transcriptf, fastaf, outfile) :
    output = {}

    with gzip.open(fastaf, 'rb') as f:
        sequence = str(f.read())

    transcript_df = pd.read_csv(transcriptf)
    transcripts = transcript_df["Transcript"].values.tolist()


    starts = transcript_df["Start"].values.tolist()
    stops = transcript_df["True_Stop"].values.tolist()

    for i in range(len(transcripts)) :
        transcript = transcripts[i]
        loc = sequence.find(transcript)
        file_name = sequence[loc:sequence.find("\\n", loc)]
        current = sequence[loc + len(file_name): sequence.find(">", loc + 1)].replace("\\n", "")

        cds = find_cds(transcript, sequence)

        if(len(cds) > 0) :
            nt = current[int(starts[i]): int(stops[i])]
            cds_nt = current[cds[0] - 1: cds[1]]

            add_dic("Sequence", nt, output)
            add_dic("AASequence", translate_sequence(nt), output)
            add_dic("CDS_AASequence", translate_sequence(cds_nt), output)
            add_dic("CDS_Start", cds[0], output)
            add_dic("CDS_Stop", cds[1], output)
        else:

            add_dic("Sequence", None, output)
            add_dic("AASequence", None, output)
            add_dic("CDS_AASequence", "", output)
            add_dic("CDS_Start", None, output)
            add_dic("CDS_Stop", None, output)

    out_file = pd.concat([pd.read_csv(transcriptf), pd.DataFrame.from_dict(output)], axis = 1)
    out_file.to_csv(outfile, index = False)

"""
Add additional analysis information to ORFs file - nucleotide sequence, 
    amino acid sequence, CDS amino acid sequence, CDS start and stop coordinates

Inputs:
- transcriptf = csv file with ORF information
- fastaf = APPRIS transcript FASTA file
- outfile = file name for output csv

Results
- outfile: csv with additional information added

"""
def main() :
    appris_file = 'appris_mouse_v2_selected.fa.gz'
    orfs = 'outputORFs.csv'
    output = 'outputORFs_comp.csv'

    write_output(orfs, appris_file, output)

main()