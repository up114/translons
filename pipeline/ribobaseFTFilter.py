
import heapq
import pandas as pd
import csv
import os
import sys
import string
import numpy as np
from ribopy import Ribo
from scipy.fft import fft, fftfreq

"""
Helper Functions
"""
###read csv file
def read_table(Experiment_file):
    table = pd.read_csv(Experiment_file)
    return table

####find studies and experiments to examine
def study_folder(database):
    study_dict = {}
    type_dict = {}

    table = read_table(database)
    study_list = table['study_name'].values.tolist()
    experiment_list = table['experiment_alias'].values.tolist()
    cell_line = table['cell_line'].values.tolist()

    for i in range(len(experiment_list)):
        if not pd.isna(study_list[i]):
            add_dic(study_list[i], experiment_list[i], study_dict)
            type_dict[experiment_list[i]] = cell_line[i]
    return study_dict, type_dict

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

####simplify transcript name to gene name
def defalternative_human_alias(x):
    x_pieces = x.split("|")
    return x_pieces[4]

####create Ribo object from path
def ribo_data(ribo_path):
    ribo_object = Ribo(ribo_path, alias=defalternative_human_alias)
    return ribo_object

####add a value to a dictionary
def add_dic(key, value, dct):
    if key in dct.keys():
        dct[key].append(value)
    else:
        dct[key] = [value]

####return p-site offset for an experiment
def psite_offset(ribo_object, exp, mmin, mmax) :
    df = (ribo_object.get_metagene("start", experiments = exp, range_lower= mmin, range_upper= mmax, sum_lengths = False, sum_references = True))

    p_site = {}

    for index, row in df.iterrows():
        max_value_index = row.iloc[35:41].idxmax() # p-site offsets are restricted to a 11-16 nt range
        offset = -1 * max_value_index

        p_site[index[1]] = offset

    return p_site

####get reads from a gene region
def get_region_reads(ribo_object, region, mmin, mmax, experiment) :
    CDS_reads = (ribo_object.get_region_counts(region_name=region,
                                          range_lower=mmin,
                                          range_upper=mmax,
                                          sum_lengths=True,
                                          sum_references=False,
                                          alias=True,
                                          experiments=experiment))
    return CDS_reads


"""
Aggregate reads across highly periodic experiments, then perform FT on ORF coordinates in aggregated dataset

Input:
- database: file location of CSV with experiments (and corresponding studies) to aggregate
- transcripts: file location of CSV with ORF coordinates

Output:
- transcript_df: contains all ORFs with FT results
- accepted_df: contains ORFs that passed the FT filter
"""
def find_FT(database, transcripts, database_path):

    transcript_df = pd.read_csv(transcripts)
    genes = transcript_df["Gene"].values.tolist()
    starts = transcript_df["Start"].values.tolist()
    stops = transcript_df["Stop"].values.tolist()

    transcript_df["Reads"] = 0
    transcript_df["Exp Reads"] = 0
    transcript_df["Max Freq"] = pd.Series([[] for _ in range(len(transcript_df))], dtype=object)
    transcript_df["FT_Accept"] = 0

    study_dict, type_dict = study_folder(database)
    exp_reads = {}
    exp_sum = 0

    # Aggregate reads
    for study in study_dict.keys():
        study_path = (database_path + "/%s" % study)  # DEPENDS ON YOUR DATABASE STRUCTURE
        if os.path.exists(study_path):
            for experiment in study_dict[study]:
                experiment_path = (study_path + "/ribo/experiments/%s.ribo" % experiment) #DEPENDS ON YOUR DATABASE STRUCTURE
                if os.path.exists(experiment_path):
                    r_file = ribo_data(experiment_path)

                    for j in r_file.experiments:  # run through each experiment
                        mmin,mmax,read_pct = intevl(experiment_path, j)
                        offset = psite_offset(r_file, j, mmin, mmax)
                        if read_pct < 0.85:
                            continue

                        exp_total = get_region_reads(r_file, "CDS", mmin, mmax, j).add(
                            get_region_reads(r_file, "UTR3_junction", mmin, mmax, j).add(
                            get_region_reads(r_file, "UTR5_junction", mmin, mmax, j).add(
                            get_region_reads(r_file, "UTR3", mmin, mmax, j).add(
                            get_region_reads(r_file, "UTR5", mmin, mmax, j)))))
                        exp_sum = exp_sum + exp_total[j].sum()


                        for k in range(mmin, (mmax + 1)):
                            df = r_file.get_coverage(experiment=experiment, range_lower=k, range_upper=k,
                                                          alias=True)
                            for l in range(len(genes)):
                                gene = genes[l]
                                start = starts[l]
                                stop = stops[l]

                                if (start <= offset[k]):
                                    continue

                                try:
                                    coverage = df[gene]
                                    offset_reads = coverage[start - offset[k]: stop - offset[k]]

                                    if (gene, start, stop) in exp_reads.keys():
                                        exp_reads[(gene,start,stop)] += offset_reads
                                    else:
                                        exp_reads[(gene,start,stop)] = offset_reads

                                except(KeyError):
                                    continue

    # Perform FT on aggregated reads
    for (gene, start, stop) in exp_reads.keys():
        signal = np.asarray(exp_reads[(gene, start, stop)], dtype=float).flatten()
        read_density = np.sum(signal) / len(signal)

        # only run FT on ORFs with at least 10% read density
        if read_density < 0.1:
            continue

        yf = abs(fft(signal))
        top_5_indices = heapq.nlargest(5, range(len(yf)), key=lambda i: yf[i])
        normalized_positions = [(i / (len(yf) - 1)) for i in top_5_indices]
        accept_flag = 1 if any(0.32 <= f <= 0.34 for f in normalized_positions) else 0

        transcript_df.loc[
            (transcript_df["Gene"] == gene) &
            (transcript_df["Start"] == start) &
            (transcript_df["Stop"] == stop),
            ["Reads", "Exp Reads", "Max Freq", "FT_Accept"]
        ] = [np.sum(signal), exp_sum, [normalized_positions], accept_flag]

    accepted_df = transcript_df[transcript_df["FT_Accept"] == 1].reset_index(drop=True)
    return transcript_df, accepted_df

"""
Writes Fourier Transform analysis results to a CSV file

Inputs:
- database: a csv file with the studies and experiments to analyze
- transcripts: a csv with gene names, start coordinate, and stop coordinate 
- database_path: RiboBase directory path

Output:
- outfile: output csv file with ORFs that pass FT filter
"""
def main():
    database = "fin_version_database_mouse_pct75.csv"
    orfs = "candidateORFs.csv"
    database_path = 'data_all'
    outfile = "outputORFs.csv"

    all_df, accepted_df = find_FT(database, orfs, database_path)
    output = pd.DataFrame.from_dict(accepted_df)
    output.to_csv(outfile, index=False)

main()