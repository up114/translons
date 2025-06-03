# translons
 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**Environment**
- ribo_env.yml: contains the necessary packages for all function other than the genomic overlap comparison
 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**Pipeline Programs**
 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**ORFPipeline.py**

This is the data pipeline used for identifying upstream and downstream open-reading frames (uORFs/dORFs) in a Ribo-Seq file. This pipeline considers each experiment independently, which was used for the EMBRYO Ribo-Seq data.
The steps of the pipeline are:
1) Identify all possible uORFs and dORFs in the APPRIS transcript annotation file
	- provide minimum and maximum bp length for uORFs
2) Determine RPF reads for all possible u/dORFs in Ribo-Seq file
3) Filter by CPM cutoff (default = 10 CPM) and reads distribution for a set of in-frame ORFs:
	- keep in-frame ORFs with a given percent of max reads in all experiments (default 	  = 75%)
	- prefer ORFs that begin with ATG
	- fallback in both cases if none meet the requirements
4) Determine 3-nucleotide periodicity via a chi-square test on global frame mapping (reads mapping to Frame 1, 2, or 3), retaining u/dORFs with p <= 0.05

Data files required
- Transcript file:
    - in Fasta GZ file format
    - Genetic sequences should be preceded with gene information formatted as: *>**ENSMUST00000070533.4**|ENSMUSG00000051951.5|OTTMUSG00000026353.2|OTTMUST00000065166.1|**Xkr4-201**|Xkr4|3634|**UTR5:1-150**|**CDS:151-2094**|UTR3:2095-3634|*
    - Sequences in CAPITAL LETTERS
    - like *appris_mouse_v2_selected.fa.gz* in references
- A .ribo ribosome profiling file with coverage data
Both the APPRIS and Ribo-Seq file can be found in *Data Files*
  

Data files outputted
- output_file: CSV with final filtered ORFs in tidy format, with all reads and periodicity outputs for each experiment
- uniqueorfs_file: CSV with the unique ORF candidates -- This is in input file for the Fourier transform

 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**ORFPipelineSummed.py**

This is the data pipeline used for identifying upstream and downstream open-reading frames (uORFs/dORFs) in a Ribo-Seq file. This pipeline sums coverage across experiments, which was used for the HIPPOCAMPAL Ribo-Seq data. 
The steps of the pipeline are:
1) Identify all possible uORFs and dORFs in the APPRIS transcript annotation file
	- provide minimum and maximum bp length for uORFs
2) Determine RPF reads for all possible u/dORFs in Ribo-Seq file
3) Filter by CPM cutoff (default = 10 CPM) and reads distribution for a set of in-frame ORFs:
	- keep in-frame ORFs with a given percent of max reads (default = 75%)
	- prefer ORFs that begin with ATG
	- fallback in both cases if none meet the requirements
4) Determine 3-nucleotide periodicity via a chi-square test on global frame mapping (reads mapping to Frame 1, 2, or 3), retaining u/dORFs with p <= 0.05

Data files required
- Transcript file:
    - in Fasta GZ file format
    - Genetic sequences should be preceded with gene information formatted as: *>**ENSMUST00000070533.4**|ENSMUSG00000051951.5|OTTMUSG00000026353.2|OTTMUST00000065166.1|**Xkr4-201**|Xkr4|3634|**UTR5:1-150**|**CDS:151-2094**|UTR3:2095-3634|*
    - Sequences in CAPITAL LETTERS
    - like *appris_mouse_v2_selected.fa.gz* in references
- A .ribo ribosome profiling file with coverage data
Both the APPRIS and Ribo-Seq file can be found in *Data Files*
  

Data files outputted
- output_file: CSV with final filtered ORFs in tidy format, with all reads and periodicity outputs -- This is in input file for the Fourier transform 


 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**ribobaseFTFilter.py**

This program aggregates reads across highly periodic experiments in RiboBase, then perform FT on ORF coordinates in aggregated dataset.


Inputs:
- database : a csv file with the studies and experiments to analyze
	- like *final_database_mouse.csv* in references, filtered by the "75periodicity_passed" column
	- these are the RiboBase studies in which >= 75% of CDS-mapped reads align to one frame
- orfs: a csv with gene names, start coordinate, and stop coordinate (column names: "Gene", "Start", "Stop")
	- this is *output_file* from ORFPipelineSummed.py and *uniqueorfs_file* from ORFPipeline.py
- database_path: RiboBase directory path

Output:
- outfile: output csv file with ORFs that pass FT filter




 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**Analysis Programs**
 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**addinfo.py**

This program adds additional information to the CSV file containing basic ORF information

Inputs:
- APPRIS file:
    - in Fasta GZ file format
    - Genetic sequences should be preceded with gene information formatted as: *>**ENSMUST00000070533.4**|ENSMUSG00000051951.5|OTTMUSG00000026353.2|OTTMUST00000065166.1|**Xkr4-201**|Xkr4|3634|**UTR5:1-150**|**CDS:151-2094**|UTR3:2095-3634|*
    - Sequences in CAPITAL LETTERS
    - like *appris_mouse_v2_selected.fa.gz* in references
- ORFs file: a CSV with transcript, start, and stop (column names: "Transcript", "Start", "True_Stop")

Output:
- outfile: csv with ORF information, with the ORF nucleotide sequence, ORF amino acid sequence, CDS amino acid sequence, and CDS start and stop coordinates added.

 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**RibobaseRNASeq.py**

This program records ORF and CDS Ribo-Seq and RNA-seq reads in RiboBase datasets

Inputs:
- database : a csv file with the studies and experiments to analyze
	- like *final_database_mouse.csv* in references
- orfs: a csv with gene names, u/dORF start and stop coordinates, and CDS start and stop coordinates
	- column names: "Gene", "Start", "Stop", "CDS_Start", "CDS_Stop"
	- this file is the output of addinfo.py, which adds the CDS coordinates
- database_path: RiboBase directory path

Output:
- outfile: output csv file with u/dORF Ribo-Seq and RNA-Seq reads in RiboBase

Note: To make a table of cell line averages for all ORFs, use the following R code:
output <- output %>% group_by(Cell_Line, Gene, Start, True_Stop, Type) %>% summarise(`ORF Reads` = mean(CPM))
output <- pivot_wider(output, names_from = Cell_Line, values_from = 'ORF Reads') %>% 
  dplyr::rename(
    Stop = True_Stop
  )
output = output 

 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
