# bTB-pipeline-comparison-data-and-tools
This is a public repository containing the processing scripts used in the bTB variant pipeline comparison. A brief description of the script is included in the first lines of the code.

Below is the process followed for the results obtained in the publication. The simulated genomes and reads have been uploaded to "zenodo.org".

1. Run "batch_reads_sim.py"
2. Input the simulated reads in each of the pipelines, adjusting the different parameters to include proximal window filtering (in the case of vSNP, BovTB or the simulated vcf files, check the directory "Proximality_Analysis") or the specific hard filter files (check the "Filter_files" directory). In order to obtain non-
parsimonious SNPs from vSNP, use the af2122_zc.vcf file included in the vSNP_processing/ directory. Note that proximal SNPs need to be removed from BovTB during output pre-processing.
4. Output pre-processing:
       a. MTBseq: go to the Amend directory and select the adequate tab file;
          [PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_[window].tab
          Use the script contained in MTBseq_vartab_processing/Amended_tab_processing/ directory and the "vcf_header.txt". You must specify the name of the 
          tab file in the code for it to recognise it. MTbseq vcf files will then be generated per sample.
       b. vSNP: go to the step2 results folder and run excel_table_processing to obtain a vcf file per position.
       c. BovTB: extract SNPs from the core genome VCF file using bcftools. If proximal SNPs are being removed, this is the moment when proximal.py should be used.
          Then use bcftools consensus caller in accordance with the parameters registered in the BovTB code and use the specific Masking files included in the work/
          directory.
5. Performance analysis: extract all of the vcf files per each of the analysis you carried out and use som.py to compare it with the simulated VCF. Note that the simulated VCF files need to be appropriately filtered as well (use bcftools, for example). Use bcftools isec to obtain FN, FP and agreeing positions.
6. Annotations: if you would like to annotate the vcf files, you may use vcf_annotation.py (requires SnpEff and the updated Mbovis annotation (https://doi.org/10.1099/acmi.0.000129). In order to annotate vcf files and obtain a summarised annotation file per sample, use annotation.py (requires sufix, for example the "0000" obtained for FN in bcftools isec).
7. FN analysis: extract positions using Allele_extraction/vcf_pos_extract.py on the FN, FP, agreeing vcf files.Then run FN_comparison.R on the csv files.
8. venn_diagram: extract all SNP positions using the script in the Allele_extraction/ directory "vcf_processing.py". Move the "[prefix]-allpositions.tsv" files to a new directory and run all_pos_ven.R to obtain the venn diagrams.
9. Tree analysis:
       a. BovTB: Use snp-sites to extract the core SNPs from the consensus genomes.
       b. Simulation: use sim_coreextractor.py on the simulated vcf files to extract the (polymorphic) core SNP alignment fasta.
       c. Run modify_samples.sh to change sample names to number codes.
       d. All: run raxml on the fasta alignment files.
       e. Homoplasies: run homoplasyFinder as indicated in its corresponding wiki.
       f. Treespace: run trees_per_pipeline.R
       g. Pairwise tree comparison: run plot_trees_wo_outgroups_final.R      
       
