
import gzip
import re
import numpy as np
import pandas as pd
from ribopy import Ribo
import os
from scipy.stats import chisquare
from statsmodels.stats.multitest import multipletests
from collections import defaultdict

"""
Helper Functions
"""

####function to create aliases and shorten transcript names
def defalternative_human_alias(x):
    x_pieces = x.split("|")
    return x_pieces[4]

####create Ribo object from path
def ribo_data(ribo_path):
    ribo_object = Ribo(ribo_path, alias=defalternative_human_alias)
    return ribo_object

####find p-site
def psite_offset(ribo_object, exp, mmin, mmax) :
    df = (ribo_object.get_metagene("start", experiments = exp, range_lower= mmin, range_upper= mmax,
                                   sum_lengths = False, sum_references = True))
    p_site = {}
    for index, row in df.iterrows() :
        max_value_index = row.iloc[35:41].idxmax()
        offset = -1 * max_value_index
        p_site[index[1]] = offset

    return p_site

#####dynamic cutoff
def intevl(experiment_path, experiment_id):
    ribo_object = ribo_data(experiment_path)
    data_tmp=ribo_object.get_length_dist("CDS")
    data_tmp.reset_index(inplace=True)
    data=data_tmp.iloc[6:26]
    pct_85=sum(data["%s"%experiment_id])*0.85
    #pct_90=sum(data["%s"%experiment_id])*0.90
    value=data[data["%s"%experiment_id]==data["%s"%experiment_id].max()]["%s"%experiment_id].values[0]
    mmin=mmax=data[data["%s"%experiment_id]==data["%s"%experiment_id].max()]['read_length'].values[0]
    while value<=pct_85 :
        if mmax<40 and mmin>21:
            if data[data['read_length']==mmax+1]["%s"%experiment_id].values[0] >= data[data['read_length']==mmin-1]["%s"%experiment_id].values[0]:
                mmax+=1
                value+=data[data['read_length']==mmax]["%s"%experiment_id].values[0]
            else:
                mmin-=1
                value+=data[data['read_length']==mmin]["%s"%experiment_id].values[0]
        elif mmax==40:
            mmin-=1
            value+=data[data['read_length']==mmin]["%s"%experiment_id].values[0]
        elif mmin==21:
            mmax+=1
            value+=data[data['read_length']==mmax]["%s"%experiment_id].values[0]
    #print(min,max)
    read_pct=value/sum(data["%s"%experiment_id]) # to check if study is suitable for analysis
    return int(mmin),int(mmax),read_pct

####Find CDS start and stop in APPRIS transcript ID
def find_CDS(transcript_name) :
    split_transcript = transcript_name.split("|")
    CDS = split_transcript[8][4:].split("-")
    return split_transcript[4], int(CDS[0]), int(CDS[1])

####add key-value pair to dictionary
def add_dic(key, value, dct):
    if key in dct.keys():
        dct[key].append(value)
    else:
        dct[key] = [value]

####find indexes of substrings in a string
def find_multiple_substrings(string, substrings):
    pattern = '|'.join(map(re.escape, substrings))
    pattern = f'(?=({pattern}))'  # Add a positive lookahead assertion
    matches = re.finditer(pattern, string)
    return [match.start(1) for match in matches]

####get reads from a gene region
def get_region_reads(ribo_object, region, mmin, mmax, experiment) :
    reads = (ribo_object.get_region_counts(region_name=region,
                                          range_lower=mmin,
                                          range_upper=mmax,
                                          sum_lengths=True,
                                          sum_references=False,
                                          alias=True,
                                          experiments=experiment))
    return reads

"""
Find all in-frame start-stop combinations

Inputs:
- min_len: min amino acid length for uORFs
- max_len: max amino acid length for uORFs
- sequence: APPRIS transcript file as a string

Results
- orf_list: contains information for all uORF/dORF candidates

Note: We only limited the length of the uORFs, not dORFs
"""
def find_ORFs(appris_file, min_len, max_len) :
    with gzip.open(appris_file, 'rt', encoding='utf-8') as f:
        sequence = (f.read())

    orf_list = {}
    file_start = sequence.find(">",0)

    # Loop through each gene in APPRIS
    while (file_start != -1):
        next_file_start = sequence.find(">", file_start + 1)
        file_name = sequence[file_start + 1 :sequence.find("\n", file_start)].strip()

        # Skip genes with no annotated 5'UTR in APPRIS transcript
        if(file_name.find('UTR5') == -1) :
            file_start = next_file_start
            continue

        # Find CDS coordinates - COULD REPLACE WITH CDS IDENTIFIED FROM ADDINFO
        gene, CDS_start, CDS_stop = find_CDS(file_name)

        current = sequence[file_start + len(file_name) + 1:next_file_start].replace("\n", "")

        # Determine indices for all start and stop codons
        starts = np.array(find_multiple_substrings(current, ["ATG", "ACG", "ATT", "CTG", "GTG", "TTG"]))
        stops = np.array(find_multiple_substrings(current, ["TAA", "TGA", "TAG"]))

        for start in starts:
            for stop in stops[stops > start] :
                stop = stop + 3
                orf_len = stop - start
                # Identify in-frame start and stop codons (length is a multiple of 3)
                if(orf_len % 3 == 0) :
                    # Identify nonoverlapping dORFs
                    if start > CDS_stop:
                        add_dic("Transcript", file_name.split("|")[0], orf_list)
                        add_dic("Gene", gene, orf_list)
                        add_dic("True_Stop", stop, orf_list)
                        add_dic("Start_Codon", current[start:start + 3], orf_list)
                        add_dic("Start", start, orf_list)
                        add_dic("Stop", stop, orf_list)
                        add_dic("Type", "nonoverlapping_dORF", orf_list)
                        break

                    # Identify uORFs: overlapping and nonoverlapping
                    if (start < CDS_start and stop < CDS_stop): # all uORFs except exclude internal ORFs
                        if (orf_len >= min_len and orf_len <= max_len) :
                            add_dic("Transcript", file_name.split("|")[0], orf_list)
                            add_dic("Gene", gene, orf_list)
                            add_dic("True_Stop", stop, orf_list)
                            add_dic("Start_Codon", current[start:start + 3], orf_list)

                            if start < CDS_start :
                                if stop < CDS_start:
                                    add_dic("Start", start, orf_list)
                                    add_dic("Stop", stop, orf_list)
                                    add_dic("Type", "nonoverlapping_uORF", orf_list)
                                else:
                                    add_dic("Start", start, orf_list)
                                    add_dic("Stop", CDS_start, orf_list)
                                    add_dic("Type", "overlapping_uORF", orf_list)
                    break
        file_start = next_file_start

    orf_df = pd.DataFrame.from_dict(orf_list)

    return orf_df


"""
Record ORF reads

Input:
- ribo_path: file location for .ribo file
- transcript_df: dataframe with the Gene, Start, Stop of the ORFs to examine
- experiments: selection of experiments in ribo file to analyze

Results:
- merged: dataframe of transcript_df with ORF reads in tidy format
"""
def orf_reads(ribo_path, transcript_df, experiments):
    genes = transcript_df["Gene"].values.tolist()
    starts = transcript_df["Start"].values.tolist()
    stops = transcript_df["Stop"].values.tolist()

    # create RiboPy objects
    r_file = Ribo(ribo_path, alias = defalternative_human_alias)

    records = []

    # for every experiment, run through the list of inputted ORFs
    for j in experiments:
        mmin, mmax, read_pct = intevl(ribo_path, j)
        offset = psite_offset(r_file, j, mmin, mmax)
        # only include experiments that pass dynamic cutoff of 85%
        if read_pct < 0.85:
            continue

        exp_total = get_region_reads(r_file, "CDS", mmin, mmax, j).add(
            get_region_reads(r_file, "UTR3_junction", mmin, mmax, j).add(
                get_region_reads(r_file, "UTR5_junction", mmin, mmax, j).add(
                    get_region_reads(r_file, "UTR3", mmin, mmax, j).add(
                        get_region_reads(r_file, "UTR5", mmin, mmax, j)))))
        exp_sum = exp_total[j].sum()

        exp_reads = defaultdict(int)
        frame1_reads = defaultdict(int)
        frame2_reads = defaultdict(int)
        frame3_reads = defaultdict(int)

        for k in range(mmin, (mmax + 1)):

            df = r_file.get_coverage(experiment=j, range_lower=k, range_upper=k,
                                     alias=True)

            for l in range(len(genes)):  # run through the gene names
                gene = genes[l]
                start = starts[l]
                stop = stops[l]

                if (start <= offset[k]):
                    continue

                try:
                    coverage = df[gene]  # store the coverage info of the gene
                    offset_reads = (coverage[start - offset[k]: stop - offset[k]])

                    exp_reads[(gene, start, stop)] += np.sum(offset_reads)
                    frame1_reads[(gene, start, stop)] += np.sum(offset_reads[::3])
                    frame2_reads[(gene, start, stop)] += np.sum(offset_reads[1::3])
                    frame3_reads[(gene, start, stop)] += np.sum(offset_reads[2::3])

                    if (k == mmax):
                        meta = transcript_df.iloc[l].to_dict()
                        record = {
                            **meta,
                            "Experiment": j,
                            "ORF Reads": exp_reads[(gene, start, stop)],
                            "Experiment Reads": exp_sum,
                            "Frames": [
                                frame1_reads[(gene, start, stop)],
                                frame2_reads[(gene, start, stop)],
                                frame3_reads[(gene, start, stop)]
                            ]
                        }
                        records.append(record)

                except(KeyError):
                    continue

    return pd.DataFrame.from_records(records)

"""
Filters uORFs by reads distribution and CPM threshold
Input:
- df: tidy DataFrame with ORF info
- experiments: list of experiment names matching df['Experiment']
- cpm_thresh: numeric CPM threshold (default 10)
- perc: fraction of group-max CPM to enforce (default 0.75)

Returns:
- tidy_out: DataFrame containing all columns that passed the filters
"""
def filter_orfs(df, experiments, cpm_thresh=10, perc=0.75):
    # copy and compute CPM
    df2 = df.copy()
    df2['CPM'] = df2['ORF Reads'] / df2['Experiment Reads'] * 1e6

    # wide summary for filtering
    summary = (
        df2
        .groupby(['Gene','Start','Stop','True_Stop','Start_Codon','Experiment'], as_index=False)
        .agg(CPM=('CPM','sum'))
        .pivot_table(
            index=['Gene','Start','Stop','True_Stop','Start_Codon'],
            columns='Experiment',
            values='CPM',
            fill_value=0
        )
        .reset_index()
    )
    cpm_cols = experiments

    # threshold any-expt CPM
    summary['max_CPM'] = summary[cpm_cols].max(axis=1)
    summary = summary[summary['max_CPM'] >= cpm_thresh]

    # containment filter - choose shortest ORF that contains all the reads in all experiments
    keep_idxs = []
    for (_, true_stop), grp in summary.groupby(['Gene','True_Stop']):
        grp_sorted = grp.sort_values('Start', ascending=False)
        kept = []
        for idx, row in grp_sorted.iterrows():
            if any((row[cpm_cols] <= grp_sorted.loc[k, cpm_cols]).all() for k in kept):
                continue
            kept.append(idx)
        keep_idxs.extend(kept)
    summary = summary.loc[keep_idxs]

    # In-frame ORF filter: ATG preference and relative abundance (CPM ≥ perc * max CPM)
    final_keys = []
    for (_, true_stop), grp in summary.groupby(['Gene','True_Stop']):
        # prefer ATG
        if (grp['Start_Codon'] == 'ATG').any():
            grp = grp[grp['Start_Codon'] == 'ATG']
        # relative abundance
        max_val = grp[cpm_cols].max()
        mask = (grp[cpm_cols] >= perc * max_val).all(axis=1)
        keep_grp = grp[mask] if mask.any() else grp
        final_keys.extend(list(zip(keep_grp['Gene'], keep_grp['Start'], keep_grp['True_Stop'])))

    key_df = pd.DataFrame(final_keys, columns=['Gene','Start','True_Stop']).drop_duplicates()
    tidy_out = df2.merge(key_df, on=['Gene','Start','True_Stop'], how='inner')

    return tidy_out.reset_index(drop=True)


"""
Output chi-square results for 3-nucleotide periodicity

Returns:
- long_df: all columns in input df with chi-square results
- filtered_periodic: columns in input df that pass chi-square
"""
def find_periodicity(df):
    long_df = df.copy()

    chi_out = find_chisquared(long_df["Frames"])

    long_df["Accept"] = chi_out["Accept"]
    long_df["Corrected p-value"] = chi_out["Corrected p-value"]

    long_df['Accept_sum'] = (
        long_df
        .groupby(['Gene', 'Start', 'Stop'])['Accept']
        .transform('sum')
    )

    # Keep only ORFs with at least one Reject==1
    filtered_periodic = long_df[long_df['Accept_sum'] > 0].copy()
    filtered_periodic.drop(columns=['Accept_sum'], inplace=True)

    return long_df, filtered_periodic


"""
Performs chi-square with FDR correction on the array containing the reads mapping 
to each ORF frame

Input:
- arr: the "Frames" column of the ORF df
Returns:
- output: dict with 
    1) "Accept" - 1 if significant, 0 if not
    2) "Corrected p-value" - FDR corrected p-value
"""
def find_chisquared(arr):
    p_vals = []
    output = {}

    for row in arr:
        total = sum(row)
        if total == 0:
            p_vals.append(np.nan)
            continue

        chi2, p_value = chisquare(row, f_exp=[total / 3, total / 3, total / 3])
        p_vals.append(p_value)

    p_vals = np.array(p_vals)
    mask = np.isfinite(p_vals)
    pval_corrected = np.empty(p_vals.shape)
    pval_corrected.fill(np.nan)

    accept_list = np.empty(p_vals.shape)
    accept_list.fill(np.nan)

    accept_list[mask], pval_corrected[mask] = multipletests(p_vals[mask], method='fdr_bh', alpha=0.05)[:2]
    output["Accept"] = accept_list.tolist()
    output["Corrected p-value"] = pval_corrected

    return output

"""
Determine ORF candidates based on experiment-dependent coverage data

Inputs:
- appris_file: APPRIS transcript annotation file in Fasta GZ file format
- ribo_file: Ribo-Seq transcript-aligned file in *.ribo format

Outputs:
- output_file: output CSV file with ORF candidate coordinates and experiment coverage data in tidy format
- uniqueorfs_file: output CVS with ORF candidate coordinates
"""
def main() :
    # Set input and output variables
    appris_file = "appris_mouse_v2_selected.fa.gz"
    ribo_file = 'merged_mapq10.ribo'
    exps = ['16cell_Control', '32cell_Control']
    output_file = 'allorfResults.csv'
    uniqueorfs_file = 'candidateORFs.csv'
    uorf_min = 33
    uorf_max = 300

    # 1) Find all ORF candidates
    allORF_df = find_ORFs(appris_file, uorf_min, uorf_max)
    # 2) Find reads for all ORF candidates
    reads_df = orf_reads(ribo_file, allORF_df, exps)
    # 3) Filter reads by reads distribution and CPM threshold
    readsFilt_df = filter_orfs(reads_df,
                               exps,
                               cpm_thresh=10,
                               perc=0.75)
    # 4) Determine 3-nucleotide periodicity with chi-square test
    full_df, sig_df = find_periodicity(readsFilt_df)


    # 5) Output final filtered ORFs in tidy format
    sig_df.to_csv(output_file, index = False)

    # 6) Extract the unique ORF candidates: remove duplicates based on key columns
    uniqueorfs_df = sig_df[
        ['Transcript', 'Gene', 'Start', 'Stop', 'Type', 'True_Stop', 'Start_Codon']].drop_duplicates()

    # 7) Output the unique ORF candidates
    uniqueorfs_df.to_csv(uniqueorfs_file, index=False)

main()