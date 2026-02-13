import pandas as pd
import csv
import os
import sys
import string
import numpy as np
from ribopy import Ribo

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
    df = read_table(database)
    study_list = df['study_name'].values.tolist()
    experiment_list = df['experiment_alias'].values.tolist()
    cell_line = df['cell_line'].values.tolist()

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
def add_dic(key, value, dict):
    if key in dict.keys():
        dict[key].append(value)
    else:
        dict[key] = [value]

####return p-site offset for an experiment
def psite_offset(ribo_object, exp, mmin, mmax) :
    df = (ribo_object.get_metagene("start", experiments = exp, range_lower= mmin, range_upper= mmax, sum_lengths = False, sum_references = True))

    p_site = {}

    for index, row in df.iterrows():
        max_value_index = row.iloc[35:41].idxmax()
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
Algorithm to record ORF and CDS Ribo-Seq and RNA-seq reads
"""
def orf_reads(database, transcripts, database_path):
    transcript_df = pd.read_csv(transcripts)

    transcript_ids = transcript_df["Transcript"].values.tolist()
    genes = transcript_df["Gene"].values.tolist()
    starts = transcript_df["Start"].values.tolist()
    stops = transcript_df["Stop"].values.tolist()

    cdsstarts = transcript_df["CDS_Start"].values.tolist()
    cdsstops = transcript_df["CDS_Stop"].values.tolist()

    study_dict, type_dict = study_folder(database)
    records = {}

    for study in study_dict.keys():
        study_path = (database_path + "/%s" % study) #DEPENDS ON YOUR DATABASE STRUCTURE
        if os.path.exists(study_path):
            for experiment in study_dict[study]:
                experiment_path = (study_path + "/ribo/experiments/%s.ribo" % experiment) #DEPENDS ON YOUR DATABASE STRUCTURE
                if os.path.exists(experiment_path):
                    r_file = ribo_data(experiment_path)

                    for j in r_file.experiments:  # run through each experiment
                        # dynamic cutoff
                        mmin,mmax,read_pct = intevl(experiment_path, j)
                        if read_pct < 0.85:
                            continue

                        # P-site offset
                        offset = psite_offset(r_file, j, mmin, mmax)

                        if r_file.has_rnaseq(experiment = j):
                            rna_seq_reads = r_file.get_rnaseq(experiments=j)
                            rna_sum = np.sum(rna_seq_reads)
                            rna_sum = np.sum(rna_sum)

                        exp_total = get_region_reads(r_file, "CDS", mmin, mmax, j).add(
                            get_region_reads(r_file, "UTR3_junction", mmin, mmax, j).add(
                            get_region_reads(r_file, "UTR5_junction", mmin, mmax, j).add(
                            get_region_reads(r_file, "UTR3", mmin, mmax, j).add(
                            get_region_reads(r_file, "UTR5", mmin, mmax, j)))))

                        exp_sum = exp_total[j].sum()

                        # only examine experiments with >100k reads
                        if exp_sum <= 100000:
                            continue

                        exp_reads = {}
                        cds_reads = {}

                        for k in range(mmin, (mmax + 1)):

                            df = r_file.get_coverage(experiment=j, range_lower=k, range_upper=k,
                                                          alias=True)

                            for l in range(len(genes)):  # run through the gene names
                                gene = genes[l]
                                start = starts[l]
                                stop = stops[l]

                                cdsstart = cdsstarts[l]
                                cdsstop = cdsstops[l]


                                if (start <= offset[k]):
                                    continue

                                try:
                                    coverage = df[gene]  # store the coverage info of the gene
                                    offset_reads = np.sum(coverage[start - offset[k]: stop - offset[k]])
                                    cds_read = np.sum(coverage[cdsstart - offset[k]: cdsstop - offset[k]])


                                    if (gene, start, stop) in exp_reads.keys():
                                        exp_reads[(gene, start, stop)] += offset_reads
                                        cds_reads[(gene, start, stop)] += cds_read
                                    else:
                                        exp_reads[(gene, start, stop)] = offset_reads
                                        cds_reads[(gene, start, stop)] = cds_read

                                    if (k == mmax):
                                        add_dic('Study', study, records)
                                        add_dic('Experiment', j, records)
                                        add_dic('Cell_Line', type_dict[j], records)
                                        add_dic('Gene', gene, records)
                                        add_dic('ORF Reads', exp_reads[(gene, start, stop)], records)
                                        add_dic('CDS Reads', cds_reads[(gene, start, stop)], records)
                                        add_dic('Experiment Reads', exp_sum, records)
                                        add_dic('Start', start, records)
                                        add_dic('Stop', stop, records)
                                        if(r_file.has_rnaseq(experiment = j)):
                                            add_dic('RNASeq_Exp', rna_sum,records)
                                            gene_rna = rna_seq_reads[rna_seq_reads.index.get_level_values(1).str.contains(transcript_ids[l])]
                                            gene_sum = np.sum(np.sum(gene_rna))
                                            add_dic('RNASeq_Gene', gene_sum, records)
                                        else:
                                            add_dic('RNASeq_Exp', None, records)
                                            add_dic('RNASeq_Gene', None, records)

                                except(KeyError):
                                    continue
    return records


"""
Determines Ribo-Seq and RNA-Seq reads for inputted ORFs

Inputs:
- database: csv file location with the study_name, experiment_alias, and cell_line of the experiments to examine
- transcripts: csv file path with the Gene, Start, Stop, CDS_Start, CDS_Stop of the ORFs to examine
- database_path: RiboBase directory path

Output:
- records: dataframe of transcript_df with ORF reads in tidy format

"""
def main() :
    database = "fin_version_database_mouse.csv"
    orfs = "outputORFs.csv"
    database_path = 'data_all'
    outfile = "outputRibobase.csv"

    output = orf_reads(database, orfs, database_path)
    output = pd.DataFrame.from_dict(output)
    output.to_csv(outfile, index=False)

main()
