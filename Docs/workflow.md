## About the page
We describe step by step analysis procedures for the manuscript about the focal amplifications in cancer. 

## Contributors

**Jake Lee**
- <https://github.com/jakelee0711>

**Lucy Jung**
- <https://github.com/YLucyJung>

## Table of contents

- [Detection of focal amplifications](#Detection-of-focal-amplifications)
- [SVs involving the boundaries of focal amplifications](#SVs-involving-the-boundaries-of-focal-amplifications)
- [Association with epigenomics data](#Association-with-epigenomics-data)
- [Association with 3D chromatin interaction data](#Association-with-3D-chromatin-interaction-data)

## Detection of focal amplifications

## SVs involving the boundaries of focal amplifications

## Association with epigenomics data
The following functions can be loaded by installing an R package TRAMP (**TR**anslocations involving oncogenic **AMP**lification) and are also in association.with.epigenomics.data.R. 

To determine which epigenomic features were associated with the initial SV events of the amplicons in breast cancers, we integrated the SVs at the amplicon boundaries with various chromatin features. This function takes the coordinates of bindings for each factor and boundary positions in 100-kb bins as inputs, computs the enrichment p-values of the factor by Fisher's exact test and displays the result.
```
association.with.chromatin.features(feature.file=featureinfile,sv.file=sv.file)
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
To analyze the association between amplicon boundaries and chromatin proximity, we obtained Hi-C data from T47D luminal breast cancer cell line which doesn't have major translocation major between the chromosomes of our interests such as between chromosomes 8, 11 and 17 and used contact frequencies normalized by balance-based method (KR normalization). This function takes amplicon boundareis and contact frequeinces in 2.5Mb as inputs and displays the association between them.
```
association.recurrence.3d.contact.t47d(rec.file=rec.file,3d.file=3d.file)
```
For the comparison of chromatin interactions between untreated- and E2-treated conditions, we first simplified translocation information connecting amplicon boundaries to an arm-level translocation network. This function takes the translocation information as an input and computed the arm-level network and visualizes.
```
translocation.network(file=file)
```
Then, the top frequently translocated chromosome pairs were determined. We obtained contact frequency information from a published, 3C-based high-throughput sequencing data in untreated- and E2-treated MCF-7 cells. The contact frequencies were combined for each chromosome arm-pair. This function takes the contact frequencies and top tranlocated arm-pairs, calculates the ratio of the arm-level contact frequencies in E2-treated cells with respect to untreated cells and compares the changes in chromatin contact with chromosome arm-level frequencies of translocations connecting the amplification boundaries. 

```
association.recurrence.3d.contact.ctr.e2.mcf7(top.trans.file=top.trans.file,3d.file=3d.file)
```
