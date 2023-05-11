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


- [Association with 3D chromatin interaction data](#Association-with-3D-chromatin-interaction-data)


## Creating metadata
Here we share the scripts for creating the R dataframes summarizing the clinical and genomic information. The three output files (`List.patients_summary.final.txt`, `List.purple.final.txt`, and `List.hrd_status.final.txt`) are available in the <a href="https://github.com/parklab/focal-amplification/blob/main/Data">Data</a> folder.

## Identifying focal amplifications
To identify focally amplified regions and to associate structural variations (SVs) to the boundaries of amplicons, we used `HMF_definition_amp_segment.R` script. This function requires the segmented allelic copy number and structural variation information produced by the HMF bioinformatic pipeline as the input as well as the dataframe summarizing the clinical information (the output from the previous section).

The output files from the `HMF_definition_amp_segment.R` were further analyzed by the scripts available in `HMF_collecting_amp_segment_descriptive_analyses.R`. Through this step, the SVs at the amplicon boundaries were summarized and annotated. Then, the copy number of the adjacent segments were analyzed to select the amplicons bordered by the unamplified segments for downsteam analyses.

Same analysis for the PCAWG dataset was performed using the `PCAWG_definition_amp_segment.R` and `PCAWG_collecting_amp_segment_descriptive_analyses`. Figure 5 was generated based on this analysis.

## Epigenomic association
The following functions can be loaded by running `association.with.epigenomics.data.R`. The relavant data files are found in the <a href="https://github.com/parklab/focal-amplification/blob/main/Data">Data</a> folder. The nessary input files for the functions are typically set as default. Figure 3 was generated based on this analysis.

To determine which epigenomic features were associated with the early SV events initiating the focal amplifications, we modelled the location of the amplicon boundaries/SVs with various epigenomic features. The function `association.with.chromatin.features` takes the coordinates of each factor and boundary positions as inputs and computes the enrichment p-values by the Fisher's exact test. 

The function `comparison.er.e2.control`compares the distributions of ERa binding intensity in E2-treated MCF7 cells and non-treated MCF7 cells and annotates major amplicon boundary hotspots.  

To study the association between the recurrence of amplicon boundaries and ERa intensity in E2-treated cells, we calculated the recurrence of patients harboring the amplification boundaries for each bin and the accumulated ERa binding intensity in the E2-treated MCF7 cells. The function `association.recurrence.e2.er.intensity` takes the information of the recurrence and the ERa intensity and displays the increase of ERa binding intensity at the binns with a high recurrence.

## HTGTS
The raw data from the HTGTS experiments is available at <ahref="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227369">GSE227369</a>. 

## Timing analysis


## Transcriptome analysis


## Data visualization
In order to visualize the patterns of SVs in the boundaries of focal amplifications, we used `Pieclust_ampseg.R` script. This accepts the summarized output txt file from the `Def_ampseg.R` script (One example result available in `Data` folder; `BreastCancer278.Ampdf.v2.full.txt`). Overall frequency of different types of boundary SVs (4 different categories including **head-to-tail SVs, fold-back inversions, translocations, and other intrachromosomal SVs**) are calculated by tumor type and the result is subject to a hierarchical clustering. The result will be visualized in pie graphs sorted by chromosome and tumor type.

To illustrate structural variations and their associated copy number information simultaneously, `SVsketch` toolkit was used (under development as an R package). Here we present an example script (`SVsketch_2-chrom_example.R`) of this toolkit to visualize complex genomic rearrangements between chromosomes 17 and 8 in patient DO1281.





## Association with 3D chromatin interaction data
The following functions can be loaded by running **association.with.epigenomics.data.R** under **focal-amplification/R** folder. The relavant data files are found in under **focal-amplification/Data** folder. The nessary input files for the functions are typically set as default.

To analyze the association between amplicon boundaries and chromatin proximity, we obtained Hi-C data from T47D luminal breast cancer cell line which doesn't have major translocation major between the chromosomes of our interests such as between chromosomes 8, 11 and 17 and used contact frequencies normalized by balance-based method (KR normalization). The function `association.recurrence.3d.contact.t47d` takes SVs from the amplicon boundareis and the location of the folder containing of contact frequeincy information in 2.5Mb as inputs and displays the association between them.

For the comparison of chromatin interactions between untreated- and E2-treated conditions in MCF7 cells, we first simplified translocation information connecting amplicon boundaries to an arm-level translocation network. The function `translocation.network` takes the translocation information as an input and computed the arm-level network and visualizes.

Then, the top frequently translocated chromosome pairs were determined. We obtained contact frequency information from a published, 3C-based high-throughput sequencing data in untreated- and E2-treated MCF-7 cells. The contact frequencies were combined for each chromosome arm-pair. The function `association.recurrence.3d.contact.ctr.e2.mcf7` takes the contact frequencies and top tranlocated arm-pairs, calculates the ratio of the arm-level contact frequencies in E2-treated cells with respect to untreated cells and compares the changes in chromatin contact with chromosome arm-level frequencies of translocations connecting the amplification boundaries. 
