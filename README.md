# awesome-deconvolution
List of deconvolution methods and resources.

# Software

## Methods

* [MuSiC](https://github.com/xuranw/MuSiC) : MuSiC is an analysis toolkit for single-cell RNA-Seq experiments

* [Bisque](https://github.com/cozygene/bisque) : An R toolkit for accurate and efficient estimation of cell composition ('decomposition') from bulk expression data with single-cell information.

* [CIBERSORT](https://cibersortx.stanford.edu/) : CIBERSORT is an analytical tool from the Alizadeh Lab and Newman Lab to impute gene expression profiles and provide an estimation of the abundances of member cell types in a mixed cell population, using gene expression data ([Newman et al 2015](https://www.nature.com/articles/nmeth.3337)).

* [CIBERSORTx](https://cibersortx.stanford.edu/) : Extension of the original CIBERSORT for RNA-seq transcriptomics data ([Newman et al 2019](https://www.nature.com/articles/s41587-019-0114-2)).

* [DeconRNASeq](http://bioconductor.org/packages/release/bioc/html/DeconRNASeq.html) : R package for deconvolution of heterogeneous tissues based on mRNA-Seq data. It modeled expression levels from heterogeneous cell populations in mRNA-Seq as the weighted average of expression from different constituting cell types and predicted cell type proportions of single expression profiles. ([Gong et al 2013](https://academic.oup.com/bioinformatics/article/29/8/1083/229442))

* [EPIC](https://epic.gfellerlab.org/) : “Estimating the Proportion of Immune and Cancer cells”. Compares the level of expression of genes in a tumor with a library of the gene expression profiles from specific cell types that can be found in tumors and uses this information to predict how many of each type of cell are present ([Racle et al 2017](https://elifesciences.org/articles/26476)).

* [MCP-counter](https://zenodo.org/record/61372#.Y-0iHXbMJPY) : Microenvironment Cell Populations-counter (MCP-counter) method, which allows the robust quantification of the absolute abundance of eight immune and two stromal cell populations in heterogeneous tissues from transcriptomic data. ([Becht et al 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1070-5))

* [ESTIMATE](https://sourceforge.net/projects/estimateproject/) :  ‘Estimation of STromal and Immune cells in MAlignant Tumours using Expression data’ (ESTIMATE)—a method that uses gene expression signatures to infer the fraction of stromal and immune cells in tumour samples. ([Yoshihara et al 2013](https://www.nature.com/articles/ncomms3612))

* [ISOpure](https://genomemedicine.biomedcentral.com/articles/10.1186/gm433) : Uses a set of tumor expression profiles and a panel of healthy tissue expression profiles to generate a purified cancer profile for each tumor sample, and an estimate of the proportion of RNA originating from cancerous cells ([Quon et al 2013](https://genomemedicine.biomedcentral.com/articles/10.1186/gm433))

* [DSA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-89#Sec4) : Digital Sorting Algorithm (DSA) for extracting cell-type specific gene expression profiles from mixed tissue samples that is unbiased and does not require prior knowledge of cell type frequencies.

## Simulation

* [SimBu](https://github.com/omnideconv/SimBu) : The goal of SimBu is to simulate pseudo-bulk RNAseq datasets with variable cell-type fractions baed on public or private single-cell RNAseq datasets.

# Resource hubs for deconvolution

* [omnideconv.org](https://omnideconv.org/) : omnideconv is an ecosystem of user-friendly tools and resources for the cell-type deconvolution of any organism and tissue profiled with bulk transcriptomics.

# See also

* [`awesome-awesomeness`](https://github.com/bayandin/awesome-awesomeness) : Large list of "awesome-*" style resources.
*  [`awesome-single-cell`](https://github.com/seandavi/awesome-single-cell) : Large list of single-cell RNA-seq resources.
