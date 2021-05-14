## About the page
We describe step by step analysis procedures for the manuscript about the focal amplifications in cancer. 

## Contributors

**Jake Lee**
- <https://github.com/jakelee0711>

**Lucy Jung**
- <https://github.com/YLucyJung>

## Table of contents

- [Defining focal amplifications and association with structural variations](#Defining-focal-amplifications-and-association-with-structural-variations)
- [Association with epigenomics data](#Association-with-epigenomics-data)
- [Association with 3D chromatin interaction data](#Association-with-3D-chromatin-interaction-data)
- [Data visualization](#Data-visualization)

## Defining focal amplifications and association with structural variations
To identify focally amplified regions and to associate structural variations (SVs) to the boundaries of amplicons, we used `Def_ampseg.R` script (current version 1.0). This function works by tumor type and requires summary information of the tumors (`Summaryinfo.table.1.19.txt`), chromosomal coordinate information (`hg19_coord.txt`), absolute copy number estimates, and annotated SV information as the input (available in `Data` folder). The output includes three key information: 1) **CNV_amplified_regions** (a segmented copy number files with annotation of focally amplified segments), 2) **SV_amp_boundaries** (a BEDPE file with annotation of types of SVs and their association with amplicons), and 3) **SV_breakpoints_long** (a similar SV file sorted by breakpoint location).

`Def_ampseg.R` first identifies baseline copy number of each chromosome arm. The baseline copy number is the most common total copy number value of the given chromosome arm. Using this baseline copy number as a reference, focally amplified region is defined in each tumor using following criteria.

- **\>3X of baseline copy number AND copy number of 6 or greater**
- **3X of baseline copy number or less but +6 copies or more from the baseline copy number**

The latter works for the focally amplified regions with amplified baseline. This happens frequently in the regions with prevalent arm-level copy gains (e.g., chromosome 8q). 

After this, the identified amplicons were merged together if the two adjacent amplicons are 1) close enough (less than 3 Mb away to each other) AND 2) the copy number of the intervening segment is clearly amplified from the baseline (2X of baseline copy number or greater AND copy number of 4 or greater). Amplicons were also filtered out if they are only moderately amplified from the adjacent segments (copy number difference less then 3). This pattern was often observed in the chromosomal regions with frequent nested or overlapped tandem duplications. Last, considering the known mechanistic relationship between the focal amplifications and chromothripsis, we expanded the amplified region boundaries if the boundary is close enough (within 1 Mb) to the copy number junction where the copy number of the adjacent segment is at the baseline copy number of given chromosome arm or less. Then SVs are associated with the boundaries of amplicons based on their physical proximity.

## Association with epigenomics data
The following functions can be loaded by running **association.with.epigenomics.data.R** under **focal-amplification/R** folder. The relavant data files are found in under **focal-amplification/Data** folder. The nessary input files for the functions are typically set as default.

To determine which epigenomic features were associated with the initial SV events of the amplicons in breast cancers, we integrated the SVs at the amplicon boundaries with various chromatin features. This function `association.with.chromatin.features` takes the coordinates of bindings for each factor and boundary positions in 100-kb bins as inputs and computs the enrichment p-values of the factor by Fisher's exact test.

The function `comparison.er.e2.control`compares the distributions of ERa binding intensity in E2-treated MCF7 cells and non-treated MCF7 cells and annotates major amplicon boundary hotspots.  

To see the association between the recurrence of amplicon boundaries and ERa intensity in E2-treated cells, we calculated the recurrence of patinets and accumulated ERa binding intensity of ERa in E2-treated MCF7 cells in 100-kb bins. This function `association.recurrence.e2.er.intensity` takes the information of the recurrence and ERa intensity and displays the increase of ERa binding intensity at the binns with a high recurrence.

## Association with 3D chromatin interaction data
The following functions can be loaded by running **association.with.epigenomics.data.R** under **focal-amplification/R** folder. The relavant data files are found in under **focal-amplification/Data** folder. The nessary input files for the functions are typically set as default.

To analyze the association between amplicon boundaries and chromatin proximity, we obtained Hi-C data from T47D luminal breast cancer cell line which doesn't have major translocation major between the chromosomes of our interests such as between chromosomes 8, 11 and 17 and used contact frequencies normalized by balance-based method (KR normalization). The function `association.recurrence.3d.contact.t47d` takes SVs from the amplicon boundareis and the location of the folder containing of contact frequeincy information in 2.5Mb as inputs and displays the association between them.

For the comparison of chromatin interactions between untreated- and E2-treated conditions in MCF7 cells, we first simplified translocation information connecting amplicon boundaries to an arm-level translocation network. The function `translocation.network` takes the translocation information as an input and computed the arm-level network and visualizes.

Then, the top frequently translocated chromosome pairs were determined. We obtained contact frequency information from a published, 3C-based high-throughput sequencing data in untreated- and E2-treated MCF-7 cells. The contact frequencies were combined for each chromosome arm-pair. The function `association.recurrence.3d.contact.ctr.e2.mcf7` takes the contact frequencies and top tranlocated arm-pairs, calculates the ratio of the arm-level contact frequencies in E2-treated cells with respect to untreated cells and compares the changes in chromatin contact with chromosome arm-level frequencies of translocations connecting the amplification boundaries. 

## Data visualization
In order to visualize the patterns of SVs in the boundaries of focal amplifications, we used `Pieclust_ampseg.R` script. This accepts the summarized output txt file from the `Def_ampseg.R` script (One example result available in `Data` folder; `BreastCancer278.Ampdf.v2.full.txt`). Overall frequency of different types of boundary SVs (4 different categories including **head-to-tail SVs, fold-back inversions, translocations, and other intrachromosomal SVs**) are calculated by tumor type and the result is subject to a hierarchical clustering. The result will be visualized in pie graphs sorted by chromosome and tumor type.

To illustrate structural variations and their associated copy number information simultaneously, `SVsketch` toolkit was used (under development as an R package). Here we present an example script (`SVsketch_2-chrom_example.R`) of this toolkit to visualize complex genomic rearrangements between chromosomes 17 and 8 in patient DO1281.
