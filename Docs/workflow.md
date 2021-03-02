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
To determine which epigenomic features were most associated with the initial SV events, we integrated the SVs at the amplicon boundaries in breast cancer samples with various chromatin features. This function takes the coordinates of bindings for each factor and boundary positions in 100 kb bins as inputs, computs the enrichment p-values of the factor by Fisher's exact test and displays the result.
```
association.with.chromatin.features<-function(file=file){...}
```
This function compares the distributions of ERa binding intensity in E2-treated MCF7 cells and non-treated MCF7 cells and annotates major amplicon boundary hotspots.  
```
comparison.er.e2.control<-function(file=file){...}
```
To see the association between the recurrence of amplicon boundaries and ERa intensity in E2-treated cells, we calculated the recurrence of patinets and accumulated ERa binding intensity of ERa in E2-treated MCF7 cells in 100-kb bins. This function takes the information of the recurrence and ERa intensity and displays the increase of ERa binding intensity at the binns with a high recurrence.
```
association.recurrence.e2.er.intensity<-function(rec.file=rec.file,er.file=er.file){...}
```
## Association with 3D chromatin interaction data
To analyze the association between amplicon boundaries and chromatin proximity, we obtained Hi-C data from T47D luminal breast cancer cell line which doesn't have major translocation major between the chromosomes of our interests such as between chromosomes 8, 11 and 17 and used contact frequencies normalized by balance-based method (KR normalization). This function takes amplicon boundareis and contact frequeinces in 2.5Mb as inputs and displays the association between them.
```
association.recurrence.3d.contact<-function(rec.file=rec.file,3d.file=3d.file){...}
```
For the comparison of chromatin interactions between untreated- and E2-treated conditions, we first simplified translocation information connecting amplicon boundaries to an arm-level translocation network. This function takes the translocation information as an input and computed the arm-level network and visualizes.
```
translocation.network<-function(file=file){...}
```


