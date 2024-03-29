# awesome-deconvolution

List of deconvolution methods and resources.

# Motivation

Deconvolution is the problem of predicting pure signals from signal mixtures. It has wide applications in many areas, from astrophysics, to image processing, to medicine, to transcriptomics. Because of this,
there are many methods out there and many more are published with each year. For anyone new or returning
to this subject, the prospect of catching up on available methods and literature can be daunting. This
repo can help with this. It can also help find connections between methods across areas, to encourage
transfer learning and reuse of existing computational tools.

# Software

## Methods

These are resources (libraries, packages, papers, etc.) which provide a specific deconvolution algorithm or approach.

* [TOAST](https://bioconductor.org/packages/release/bioc/html/TOAST.html) : TOols for the Analysis of heterogeneouS Tissues. Tools for the analysis of heterogeneous tissues. (see also: [Li and Wu 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1778-0)); tags: bioconductor; transcriptomics; reference-free; partial_reference-free

* [TOAST-csDeconv])(https://bioconductor.org/packages/release/bioc/html/TOAST.html) : Function to improve the feature selection in reference-free deconvolution through cross-cell type differential analysis. (see also: [Li and Wu 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1778-0)) tags: methylation; dna_methylation; dnam; epigenetics

* [SingleR](https://bioconductor.org/packages/release/bioc/html/SingleR.html) : Performs unbiased cell type recognition from single-cell RNA sequencing data, by leveraging reference transcriptomic datasets of pure cell types to infer the cell of origin of each single cell independently (see also: [Aran et al 2019](https://www.nature.com/articles/s41590-018-0276-y)); tags: single-cell; transcriptomics; bioconductor; r

* [xCell](https://xcell.ucsf.edu/) : xCell is a webtool that performs cell type enrichment analysis from gene expression data for 64 immune and stroma cell types (see also: [Aran et al 2017](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1349-1)); tags: bulk; transcriptomics

* [DCQ](http://dcq.tau.ac.il/) : Digital cell quantification (DCQ). Server for digital cell quantity estimation (see also: [Altboum et al 2014](https://www.embopress.org/doi/full/10.1002/msb.134947)); tags: immune; immune_cell; murine; mouse;

* [seq-ImmuCC](http://218.4.234.74:3200/immune/) : ImmuCC server was developed to estimate the relative immune cell compositions in mouse tissue from the transcriptomal data profiled on both microarray and RNA-Seq platforms (see also: (Chen et al 2018)[https://www.frontiersin.org/articles/10.3389/fimmu.2018.01286/full]); tags: immune; immune_cell; murine; mouse; transcriptomics

* [mMCP-counter](https://github.com/cit-bioinfo/mMCP-counter) : Murine version of MCP-counter, a tool to estimate the immune and stromal composition of heterogeneous tissue, from transcriptomic data. It is distributed as a R package. (see also: (Petitprez et al 2020)[https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00783-w]); tags: murine; mouse; r; mcp

* [MCPcounter](https://github.com/ebecht/MCPcounter) :  The Microenvironment Cell Populations-counter (MCP-counter) method, which allows the robust quantification of the absolute abundance of eight immune and two stromal cell populations in heterogeneous tissues from transcriptomic data (see also: [Becht et al 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1070-5)); tags: immune; immune_cell; blood; blood_cell; transcriptomics; microenvironment; human

* [TIMER](http://cistrome.org/TIMER/) : Tumor IMmune Estimation Resource. TIMER is a web resource for systematical evaluations of the clinical impact of different immune cells in diverse cancer types; tags: immune; immune_system; blood; blood_cell; tumor; cancer

* [quanTIseq](https://icbi.i-med.ac.at/software/quantiseq/doc/index.html) : quanTIseq is a computational pipeline for the quantification of the Tumor Immune contexture from human RNA-seq data (see also: [Finotello et al 2019](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0638-6)); tags: immune; immune_system; tumor; cancer; 

* [SPOTlight](https://marcelosua.github.io/SPOTlight/) : SPOTlight provides a tool that enables the deconvolution of mixtures of cells from a single-cell reference (see also: [Bayes et al 2021](https://marcelosua.github.io/SPOTlight/)); tags: r; bioconductor; spatial; transcriptomics; spot; spot_deconvolution 

* [Tangram](https://github.com/broadinstitute/Tangram) : Tangram is a Python package, written in PyTorch and based on scanpy, for mapping single-cell (or single-nucleus) gene expression data onto spatial gene expression data. (see also: [Biancalani and Scalia et al 2021](https://www.nature.com/articles/s41592-021-01264-7)); tags: python; pytorch; scanpy; spatial; transcriptiomics; spot; spot_deconvolution

* [MuSiC](https://github.com/xuranw/MuSiC) : MUlti-Subject SIngle Cell deconvolution (MuSiC) is an analysis toolkit for single-cell RNA-Seq experiments ([Wang et al 2019](https://www.nature.com/articles/s41467-018-08023-x); tags: single-cell; transcriptomics; rna-seq; bulk).

* [MuSiC2](https://github.com/Jiaxin-Fan/MuSiC2) : MUlti-Subject SIngle Cell deconvolution 2 (MuSiC2). Cell type deconvolution for multi-condition bulk RNA-seq data ([Fan et al 2022](https://academic.oup.com/bib/article-abstract/23/6/bbac430/6751147?redirectedFrom=fulltext); tags: single-cell; transcriptomics; rna-seq; bulk; case_control)

* [Bisque](https://github.com/cozygene/bisque) : An R toolkit for accurate and efficient estimation of cell composition ('decomposition') from bulk expression data with single-cell information ([Jew et al 2020](https://www.nature.com/articles/s41467-020-15816-6); tags: single-cell; transcriptomics; rna-seq; bulk).

* [CIBERSORT](https://cibersortx.stanford.edu/) : CIBERSORT is an analytical tool from the Alizadeh Lab and Newman Lab to impute gene expression profiles and provide an estimation of the abundances of member cell types in a mixed cell population, using gene expression data ([Newman et al 2015](https://www.nature.com/articles/nmeth.3337)).

* [CIBERSORTx](https://cibersortx.stanford.edu/) : Extension of the original CIBERSORT for RNA-seq transcriptomics data ([Newman et al 2019](https://www.nature.com/articles/s41587-019-0114-2)).

* [EpiDISH](https://github.com/sjczheng/EpiDISH) : Epigenetic Dissection of Intra-Sample-Heterogeneity. (see also: [Teschendorff et al 2017](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1511-5)); tags: methylation; dna_methylation; dnam; epigenetics;

* [MethylCIBERSORT](https://zenodo.org/record/1298968) : Reference-based method using CBS algorithm. Implemented as an R package (see also: [Chakravarthy et al 2018](https://www.nature.com/articles/s41467-018-05570-1#Sec12)); tags: methylation; dna_methylation; dnam; epigenetics;

* [MethylResolver](https://github.com/darneson/MethylResolver) : Robust method for deconvolving bulk tissue methylation data using least trimmed squares (LTS) regression (see also: [](https://www.nature.com/articles/s42003-020-01146-2)); tags: methylation; dna_methylation; dnam; epigenetics; r

* [ARIC](https://xwanglabthu.github.io/ARIC/) : RAccurate and robust inference of cell type proportions from bulk gene expression or DNA methylation data (see also: [Zhang et al 2022](https://academic-oup-com.proxy1.library.jhu.edu/bib/article/23/1/bbab362/6361035)). tags: transcriptomics; methylation; dna_methylation; dnam; epigenetics; python; svr;

* [RefFreeEWAS](https://cran.r-project.org/src/contrib/Archive/RefFreeEWAS/) : Reference-free method for conducting EWAS while deconvoluting DNA methylation arising as mixtures of cell types.  This method is similar to surrogate variable analysis (SVA and ISVA), except that it makes additional use of a biological mixture assumption. (see also: [Houseman et al 2016](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1140-4#Sec10)). tags: methylation; dna_methylation; dnam; epigenetics; r; cran; toast; 

* [Houseman et al 2012](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-86) : Early and influential reference-based method for deconvolutino of blood cell types from DNA methylation array data. (see also: [minfi](https://www.bioconductor.org/packages/release/bioc/html/minfi.html)). tags: methylation; dna_methylation; dnam; epigenetics; minfi; r; bioconductor

* [DeconRNASeq](http://bioconductor.org/packages/release/bioc/html/DeconRNASeq.html) : R package for deconvolution of heterogeneous tissues based on mRNA-Seq data. It modeled expression levels from heterogeneous cell populations in mRNA-Seq as the weighted average of expression from different constituting cell types and predicted cell type proportions of single expression profiles. ([Gong et al 2013](https://academic.oup.com/bioinformatics/article/29/8/1083/229442))

* [EPIC](https://epic.gfellerlab.org/) : “Estimating the Proportion of Immune and Cancer cells”. Compares the level of expression of genes in a tumor with a library of the gene expression profiles from specific cell types that can be found in tumors and uses this information to predict how many of each type of cell are present ([Racle et al 2017](https://elifesciences.org/articles/26476)).

* [MCP-counter](https://zenodo.org/record/61372#.Y-0iHXbMJPY) : Microenvironment Cell Populations-counter (MCP-counter) method, which allows the robust quantification of the absolute abundance of eight immune and two stromal cell populations in heterogeneous tissues from transcriptomic data. ([Becht et al 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1070-5))

* [ESTIMATE](https://sourceforge.net/projects/estimateproject/) :  ‘Estimation of STromal and Immune cells in MAlignant Tumours using Expression data’ (ESTIMATE)—a method that uses gene expression signatures to infer the fraction of stromal and immune cells in tumour samples. ([Yoshihara et al 2013](https://www.nature.com/articles/ncomms3612))

* [ISOpure](https://genomemedicine.biomedcentral.com/articles/10.1186/gm433) : Uses a set of tumor expression profiles and a panel of healthy tissue expression profiles to generate a purified cancer profile for each tumor sample, and an estimate of the proportion of RNA originating from cancerous cells ([Quon et al 2013](https://genomemedicine.biomedcentral.com/articles/10.1186/gm433))

* [DSA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-89#Sec4) : Digital Sorting Algorithm (DSA) for extracting cell-type specific gene expression profiles from mixed tissue samples that is unbiased and does not require prior knowledge of cell type frequencies.

* [dtangle](https://cran.r-project.org/web/packages/dtangle/index.html) : Deconvolving cell types from high-throughput gene profiling data ([Hunt et al 2019](https://academic.oup.com/bioinformatics/article/35/12/2093/5165376?login=false)).

* [Bulk2Space](https://github.com/ZJUFanLab/bulk2space) : A spatial deconvolution method based on deep learning frameworks, which converts bulk transcriptomes into spatially resolved single-cell expression profiles ([Liao et al 2022](https://www.biorxiv.org/content/10.1101/2022.01.15.476472v1)).

* [cell2location](https://github.com/BayraktarLab/cell2location/) : Comprehensive mapping of tissue cell architecture via integrated single cell and spatial transcriptomics ([Kleshchevnikov et al 2022](https://www.nature.com/articles/s41587-021-01139-4)).

* [BayesSpace](http://www.bioconductor.org/packages/release/bioc/html/BayesSpace.html) : Clustering and Resolution Enhancement of Spatial Transcriptomes ([Zhao et al 2021](https://www.nature.com/articles/s41587-021-00935-2)).

* [ImmuCC](https://github.com/wuaipinglab/ImmuCC) : Cell-Centric View of Tissue Transcriptome Measuring Cellular Compositions of Immune Microenvironment From Mouse RNA-Seq Data ([Chen et al 2017](https://www.nature.com/articles/srep40508); tags: immune cells; blood cells; mouse; rna-seq).

* [SCDC](https://meichendong.github.io/SCDC/) : Bulk Gene Expression Deconvolution by Multiple Single-Cell RNA Sequencing References ([Dong et al 2021](https://academic.oup.com/bib/article/22/1/416/5699815); tags: single-cell; rna-seq; transcriptomics).

* [SPLITR](https://www.biorxiv.org/content/10.1101/2021.01.21.426000v1) : Integrates single-nucleus and bulk RNA-seq data, enabling phenotype-aware deconvolution and correcting for systematic discrepancies between bulk and single-cell data ([Park et al 2021](https://www.biorxiv.org/content/10.1101/2021.01.21.426000v1); tags: eqtl; gwas; genetics).

# Frameworks, workflows, simulations, and benchmarks

These are resouces which provide a variety of different tools and utlities for conducting a deconvolution experiment.
These may faciliate handling of simulations, generation of pseudobulks, harmonization of multiple references or
matched and mismatched bulk and single-cell data, quality-control, and more.

* [immunedeconv](https://github.com/omnideconv/immunedeconv) : An R package for unified access to computational methods for estimating immune cell fractions from bulk RNA sequencing data (see also: [Sturm et al 2019](https://academic.oup.com/bioinformatics/article/35/14/i436/5529146)); tags: R; Rstats; Rpackage; immune; immune_system; blood; blood_cells; framework

* [lute](https://github.com/metamaden/lute) : An R package for bulk deconvolution access, simulation, optimization, and benchmarking. tags: R; Rstats; bioconductor; framework; size_factors; music; bisque; deconrnaseq; epic; nnls

* [Hippen et al 2023](https://github.com/greenelab/deconvolution_pilot/tree/main/scripts/deconvolution) : Snakemake workflow for deconvolution of bulk tumor tissues (see also: [Hippen et al 2023](https://www.biorxiv.org/content/10.1101/2022.12.04.519045v2)). tags: snakemake; workflow; pipeline; tumor; cancer; bulk

* [SimBu](https://github.com/omnideconv/SimBu) : The goal of SimBu is to simulate pseudo-bulk RNAseq datasets with variable cell-type fractions baed on public or private single-cell RNAseq datasets ([Dietrich et al 2022](https://academic.oup.com/bioinformatics/article/38/Supplement_2/ii141/6702009); tags: pseudobulk; single-cell; rna-seq; simulation; bias; cell_size; scale_factor; transcriptomics; bulk).

* [splatter](http://bioconductor.org/packages/release/bioc/html/splatter.html) : Splatter is a package for the simulation of single-cell RNA sequencing count data ([Zappia et al 2017](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0); tags: simulation; single-cell; rna-seq; bias; transcriptomics).

* [TME](https://olliecast.shinyapps.io/Deconvolution_Benchmarking/) : Benchmarking of Tumour Microenvrionment Cell Type Estimation From Bulk RNA. (see also: [Jiménez-Sánchez et al 2018](https://www.biorxiv.org/content/10.1101/437533v2.full)); tags: benchmark; tumor; cancer; microenvironment

* [SCDC](https://meichendong.github.io/SCDC/) : Bulk Gene Expression Deconvolution by Multiple Single-Cell RNA Sequencing References ([Dong et al 2021](https://academic.oup.com/bib/article/22/1/416/5699815); tags: single-cell; rna-seq; transcriptomics).

* [TOAST](https://bioconductor.org/packages/release/bioc/html/TOAST.html) : TOols for the Analysis of heterogeneouS Tissues. Tools for the analysis of heterogeneous tissues. (see also: [Li and Wu 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1778-0)); tags: bioconductor; transcriptomics; reference-free; partial_reference-free

# Resource hubs for deconvolution

These are cetralized repositories for well-annotated datasets and other entities relevant to deconvolution and 
benchmarking of deconvolution algorithms.

* [omnideconv.org](https://omnideconv.org/) : omnideconv is an ecosystem of user-friendly tools and resources for the cell-type deconvolution of any organism and tissue profiled with bulk transcriptomics.

# See also

If you found `awesome-deconvolution` helpful, you may also be interested in the following:

* [`awesome-awesomeness`](https://github.com/bayandin/awesome-awesomeness) : Large list of "awesome-*" style resources.
*  [`awesome-single-cell`](https://github.com/seandavi/awesome-single-cell) : Large list of single-cell RNA-seq resources.
