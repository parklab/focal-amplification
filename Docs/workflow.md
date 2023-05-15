## Contributors

**Jake Lee**
- <https://github.com/jakelee0711>

**Lucy Jung**
- <https://github.com/YLucyJung>

## Table of contents

- [01_Creating_metadata](#Creating-metadata)
- [02_Identifying_focal_amplifications](#Identifying-focal-amplifications)
- [03_Epigenomic_association](#Epigenomic-association)
- [04_HTGTS](#HTGTS)
- [05_Timing_analysis](#Timing-analysis)
- [06_Transcriptome_analysis](#Transcriptome-analysis)
- [07_Data_visualization](#Data-visualization)

The computational scripts described here are available in the [R](https://github.com/parklab/focal-amplification/blob/main/R) folder.

## Creating metadata
Here we share the scripts for creating the R dataframes summarizing the clinical and genomic information. The three output files (`List.patients_summary.final.txt`, `List.purple.final.txt`, and `List.hrd_status.final.txt`) are available in the <a href="https://github.com/parklab/focal-amplification/blob/main/Data">Data</a> folder.

## Identifying focal amplifications
To identify focally amplified regions and to associate structural variations (SVs) to the boundaries of amplicons, we used `HMF_definition_amp_segment.R` script. This function requires the segmented allelic copy number and structural variation information produced by the HMF bioinformatic pipeline as the input as well as the dataframe summarizing the clinical information (the output from the previous section).

The output files from the `HMF_definition_amp_segment.R` were further analyzed by the scripts available in `HMF_collecting_amp_segment_descriptive_analyses.R`. Through this step, the SVs at the amplicon boundaries were summarized and annotated. Then, the copy number of the adjacent segments were analyzed to select the amplicons bordered by the unamplified segments for downsteam analyses.

Same analysis for the PCAWG dataset was performed using the `PCAWG_definition_amp_segment.R` and `PCAWG_collecting_amp_segment_descriptive_analyses`. Figure 5 was generated based on this analysis.

## Epigenomic association
The following functions can be loaded by running `association.with.epigenomics.data.R`. The relavant data files are found in the <a href="https://github.com/parklab/focal-amplification/blob/main/Data">Data</a> folder. The nessary input files for the functions are typically set as default. Figure 3 was generated based on this analysis.

To determine which epigenomic features were associated with the early SV events initiating the focal amplifications, we modelled the location of the amplicon boundaries/SVs with various epigenomic features. The function `association.with.chromatin.features` takes the coordinates of each factor and boundary positions as inputs and computes the enrichment p-values by the Lasso regression. 

The function `comparison.er.e2.control`compares the distributions of ERa binding intensity in E2-treated MCF7 cells and non-treated MCF7 cells and annotates major amplicon boundary hotspots.  

To study the association between the recurrence of SV breakpoints and ERa binding in E2-treated cells in unamplified regions, we calculated the recurrence of patients harboring SV breakpoints for each bin and the ERa binding in the E2-treated MCF7 cells. The function `association.recurrence.e2.er.non.amp` takes the information of the recurrencee, the ERa binding, and displays the increase of the percentage of regions with ERa binding with higher recurrences.

To analyze the association between break points by ER treatment by HTGTS experiments (see below) and epigenomic features, 
we modeled the location of the breakpoints with various epigenomic features. The function `association.with.chromatin.features.htgts` takes the coordinates of each factor and HTGTS breakpoint positions as inputs and computes the enrichment p-values by the Lasso regression. 

## HTGTS
The raw data from the HTGTS experiments is available at [GEO GSE227369](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227369). The raw data was primarily analyzed using the scripts in the `HTGTS_data_processing.R`. The main output file from this script, `HTGTS.lograt.average.allcells.twotarget.gene.v2.corrected.allinfo.txt`, is available in the <a href="https://github.com/parklab/focal-amplification/blob/main/Data">Data</a> folder.

The scripts used for downstream analyses of the HTGTS dataset are available in the `HTGTS_annotation_visualization.R`. The raw output from the GSEA analysis, `GSEA.report_for_na_pos_1650575918912.tsv`, is available in the <a href="https://github.com/parklab/focal-amplification/blob/main/Data">Data</a>. Relevant for Figure 3.

## Timing analysis
To analyze the timing of copy-number gains in the breast cancer genomes, we used two different methods, relative and absolute timing of the segments as described in the manuscript. First, relative timing was analyzed based on the methods used in the PCAWG analysis (`MutationTimeR` package was used) using our custom code `TimeR_CNA_WGD.R`. Synchronicity of the copy-number gains were assessed using `TimeR_synchronicity_analysis.R`.

For absolute timing, we determined the burden of pre-amplification mutation for each amplified genomic segment using `VEClonal_calculation_binomial.R`. Using the output of this analysis, we estimated the timing of major copy-number gains as well as the non-bridge arm gains from the select cases with TB amplification based on the scripts collected in `VEClonal_timing_analysis.R`. Figure 4 was generated based on this analysis.

## Transcriptome analysis
We devised a score reflecting mRNA expression of the estrogen-responsive genes and applied this to the 263 breast cancer cases with available RNA sequencing datasets. The scripts used in this analysis are available in `ERalpha_transcriptome_activity.R`. Relevant for Figure 4.

We also analyzed the expression of amplified genes and correlated with their knockout phenotype in the CRISPR screen data. Related scripts are available in `Amplified_genes_expression.R`.

## Data visualization
To illustrate structural variations and their associated copy number information, We used `SVsketch_20chrom_HMF_clean.R`. This was used in Figures 1, 2, and 4.

To visualize the genomic rearrangement landscape in 780 breast cancers, we used `Plot_SV_matrix.R` script. Relevant for Figure 1.

To visualize the amplified regions, their boundaries, and the association with the HTGTS breakpoints, we used `Plot_amplicon_HTGTS.R` script.

We used `Oncoprint_clinical_signatures.R` script to illustrate the genomic alteration landscape (Extended Data Fig. 2d). The output files were integrated using the Adobe Illustrator.




