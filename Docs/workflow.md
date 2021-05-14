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
To identify focally amplified regions and to associate structural variations (SVs) to the boundaries of amplicons, we used `Def_ampseg.R` script (current version 1.0). This function works by tumor type and requires summary information of the tumors (Summaryinfo.table.1.19.txt), chromosomal coordinate information (hg19_coord.txt), absolute copy number estimates, and annotated SV information as the input (available in Data folder). The output includes three key information: 1) **CNV_amplified_regions** (a segmented copy number files with annotation of focally amplified segments), 2) **SV_amp_boundaries** (a BEDPE file with annotation of types of SVs and their association with amplicons), and 3) **SV_breakpoints_long** (a similar SV file sorted by breakpoint location).

`Def_ampseg.R` first identifies baseline copy number of each chromosome arm. The baseline copy number is the most common total copy number value of the given chromosome arm. Using this baseline copy number as a reference, focally amplified region is defined in each tumor using following criteria.

- **\>3X of baseline copy number AND copy number of 6 or greater**
- **3X of baseline copy number or less but +6 copies or more from the baseline copy number**

The latter works for the focally amplified regions with amplified baseline. This happens frequently in the regions with prevalent arm-level copy gains (e.g., chromosome 8q). 

After this, the identified amplicons were merged together if the two adjacent amplicons are 1) close enough (less than 3 Mb away to each other) AND 2) the copy number of the intervening segment is clearly amplified from the baseline (2X of baseline copy number or greater AND copy number of 4 or greater). Amplicons were also filtered out if they are only moderately amplified from the adjacent segments (copy number difference less then 3). This pattern was often observed in the chromosomal regions with frequent nested or overlapped tandem duplications. Last, considering the known mechanistic relationship between the focal amplifications and chromothripsis, we expanded the amplified region boundaries if the boundary is close enough (within 1 Mb) to the copy number junction where the copy number of the adjacent segment is at the baseline copy number of given chromosome arm or less. Then SVs are associated with the boundaries of amplicons based on their physical proximity.

## Association with epigenomics data
The following functions can be loaded by installing an `R package TRAMP` (**TR**anslocations involving oncogenic **AMP**lification) and are also in `association.with.epigenomics.data.R` under `focal-amplification/R` folder. The relavant data files are found in under `focal-amplification/Data` folder.

To determine which epigenomic features were associated with the initial SV events of the amplicons in breast cancers, we integrated the SVs at the amplicon boundaries with various chromatin features. This function takes the coordinates of bindings for each factor and boundary positions in 100-kb bins as inputs and computs the enrichment p-values of the factor by Fisher's exact test.
```
association.with.chromatin.features(feature.file=feature.file,sv.file=sv.file)
```
This function compares the distributions of ERa binding intensity in E2-treated MCF7 cells and non-treated MCF7 cells and annotates major amplicon boundary hotspots.  
```
comparison.er.e2.control(file=file)
```
To see the association between the recurrence of amplicon boundaries and ERa intensity in E2-treated cells, we calculated the recurrence of patinets and accumulated ERa binding intensity of ERa in E2-treated MCF7 cells in 100-kb bins. This function takes the information of the recurrence and ERa intensity and displays the increase of ERa binding intensity at the binns with a high recurrence.
```
association.recurrence.e2.er.intensity(rec.file=rec.file,er.file=er.file)
```
## Association with 3D chromatin interaction data
To analyze the association between amplicon boundaries and chromatin proximity, we obtained Hi-C data from T47D luminal breast cancer cell line which doesn't have major translocation major between the chromosomes of our interests such as between chromosomes 8, 11 and 17 and used contact frequencies normalized by balance-based method (KR normalization). This function takes SVs from the amplicon boundareis and the location of the folder containing of contact frequeincy information in 2.5Mb as inputs and displays the association between them.
```
association.recurrence.3d.contact.t47d(sv.file=sv.file,3d.file.folder=3d.file.folder)
```
For the comparison of chromatin interactions between untreated- and E2-treated conditions, we first simplified translocation information connecting amplicon boundaries to an arm-level translocation network. This function takes the translocation information as an input and computed the arm-level network and visualizes.
```
translocation.network(file=file)
```
Then, the top frequently translocated chromosome pairs were determined. We obtained contact frequency information from a published, 3C-based high-throughput sequencing data in untreated- and E2-treated MCF-7 cells. The contact frequencies were combined for each chromosome arm-pair. This function takes the contact frequencies and top tranlocated arm-pairs, calculates the ratio of the arm-level contact frequencies in E2-treated cells with respect to untreated cells and compares the changes in chromatin contact with chromosome arm-level frequencies of translocations connecting the amplification boundaries. 

```
association.recurrence.3d.contact.ctr.e2.mcf7(top.trans.file=top.trans.file,3d.file=3d.file)
```
## Data visualization
