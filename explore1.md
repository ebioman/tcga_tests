TCGA tests
================

``` r
library(devtools)
```

    ## Loading required package: usethis

``` r
# get clinical data 
library(TCGAbiolinks)
library(janitor)
```

    ## 
    ## Attaching package: 'janitor'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     chisq.test, fisher.test

``` r
library(SummarizedExperiment)
```

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following objects are masked from 'package:GenomicRanges':
    ## 
    ##     intersect, setdiff, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following object is masked from 'package:matrixStats':
    ## 
    ##     count

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(here)
```

    ## here() starts at /workspaces/TCGA

``` r
library(rmarkdown)
library(org.Hs.eg.db)
```

    ## Loading required package: AnnotationDbi

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## 

``` r
library("clusterProfiler")
```

    ## 

    ## clusterProfiler v4.14.6 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
    ## 
    ## Please cite:
    ## 
    ## S Xu, E Hu, Y Cai, Z Xie, X Luo, L Zhan, W Tang, Q Wang, B Liu, R Wang,
    ## W Xie, T Wu, L Xie, G Yu. Using clusterProfiler to characterize
    ## multiomics data. Nature Protocols. 2024, 19(11):3292-3320

    ## 
    ## Attaching package: 'clusterProfiler'

    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(ggplot2)
library(ggfortify)
library(ComplexHeatmap)
```

    ## Loading required package: grid

    ## ========================================
    ## ComplexHeatmap version 2.22.0
    ## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
    ## Github page: https://github.com/jokergoo/ComplexHeatmap
    ## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
    ## 
    ## If you use it in published research, please cite either one:
    ## - Gu, Z. Complex Heatmap Visualization. iMeta 2022.
    ## - Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
    ##     genomic data. Bioinformatics 2016.
    ## 
    ## 
    ## The new InteractiveComplexHeatmap package can directly export static 
    ## complex heatmaps into an interactive Shiny app with zero effort. Have a try!
    ## 
    ## This message can be suppressed by:
    ##   suppressPackageStartupMessages(library(ComplexHeatmap))
    ## ========================================

``` r
library(maftools)
library(httpgd)
library(survival)
```

# PCA of TCGA expression data

This follows pretty much the blog post from
[Chatomics1](https://divingintogeneticsandgenomics.com/post/pca-tcga/)
and
[Chatomics2](https://divingintogeneticsandgenomics.com/post/pca-tcga2/?s=09)
since I wanted to test a few things and that was going right into a
similar direction :) The gist is essentially how to access TCGA through
R api and doing exploratory data anaylsis (EDA). For this purpose, we
download data from 2 different cancer types, Lung adenocarcinoma (LUAD)
and LUSC (lung squamous cell carcinoma).

There is as well from TCGA a nice
[tutorial](https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/analysis.html)
which we might want to use to cycle back or re-investigate again
occasionally.

Similarly, there is [here a case
study](https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/casestudy.html)
for pan-cancer BRCA downstram analysis which might be very helpful.

## LUAD

Get Lung cancer data , more detailed information for downloading can be
found [in the
doscs](https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html)

``` r
# This one should only be done once and the rds file been kept
# Therefore it is excluded from the report generation
LUAD_query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")

GDCdownload(LUAD_query)
TCGA_LUAD_data <- GDCprepare(LUAD_query)
saveRDS(TCGA_LUAD_data, "./data/TCGA_LUAD_SummarizedExperiment.rds")
```

``` r
TCGA_LUAD_data<- readRDS("./data/TCGA_LUAD_SummarizedExperiment.rds")
# save the raw counts matrix 
TCGA_LUAD_mat<- assay(TCGA_LUAD_data)

methy_LUAD <- colData(TCGA_LUAD_data) %>% as.data.frame() %>% tabyl(`paper_CIMP.methylation.signature.`)
```

We have in total 600 samples in the sunmmarized experiment Lets have a
look at methylation sub-types

We can see that there is plenty of NA values. In the raw count data,
they have ensemble IDs for each gene :

## LUSC

Get Lung cancer data , more detailed information for downloading can be
found [in the
doscs](https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html)

``` r
# This one should only be done once and the rds file been kept
# Therefore it is excluded from the report generation
LUSC_query <- GDCquery(project = "TCGA-LUSC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")

GDCdownload(LUSC_query)
TCGA_LUSC_data <- GDCprepare(LUSC_query)
saveRDS(TCGA_LUSC_data, "./data/TCGA_LUSC_SummarizedExperiment.rds")
```

``` r
TCGA_LUSC_data<- readRDS("./data/TCGA_LUSC_SummarizedExperiment.rds")
# save the raw counts matrix 
TCGA_LUSC_mat<- assay(TCGA_LUSC_data)
```

We have in total 562 samples in the sunmmarized experiment

In the raw count data, they have ensemble IDs for each gene :

## Harmonization

We want to have instead of ensemble IDs the GENE ID or better even gene
symbols.

``` r
TCGA_LUAD_genes <- rownames(TCGA_LUAD_mat)%>% tibble::enframe() %>% mutate(ENSEMBL=stringr::str_replace(value,pattern="\\.[0-9]+",replace=""))
TCGA_LUSC_genes <- rownames(TCGA_LUSC_mat)%>% tibble::enframe() %>% mutate(ENSEMBL=stringr::str_replace(value,pattern="\\.[0-9]+",replace=""))
```

These have (non-surprisingly) the same sizes of total genes: 60660

Now there are some which are duplicated in there on the gene level as we
have multiple transcripts per genes. They used a package called
[clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
for the conversion of the IDs, I think personally I would have rather
gone for
[biomart](https://bioconductor.org/packages/release/bioc/html/biomaRt.html).
I assume because this package allows more complex analysis on top of it
and therefore reduces the number of total dependencies [see recent
paper](https://www.cell.com/the-innovation/fulltext/S2666-6758(21)00066-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2666675821000667%3Fshowall%3Dtrue).

``` r
example_dup <- clusterProfiler::bitr(TCGA_LUAD_genes$ENSEMBL, 
                      fromType = "ENSEMBL",
                      toType = "SYMBOL",
                      OrgDb = org.Hs.eg.db) %>%
        janitor::get_dupes(SYMBOL) %>%
        head()
```

    ## 'select()' returned 1:many mapping between keys and columns

    ## Warning in clusterProfiler::bitr(TCGA_LUAD_genes$ENSEMBL, fromType = "ENSEMBL",
    ## : 40.14% of input gene IDs are fail to map...

``` r
TCGA_LUAD_gene_map <- clusterProfiler::bitr(TCGA_LUAD_genes$ENSEMBL, 
                      fromType = "ENSEMBL",
                      toType = "SYMBOL",
                      OrgDb = org.Hs.eg.db) %>%
                distinct(SYMBOL, .keep_all = TRUE)
```

    ## 'select()' returned 1:many mapping between keys and columns

    ## Warning in clusterProfiler::bitr(TCGA_LUAD_genes$ENSEMBL, fromType = "ENSEMBL",
    ## : 40.14% of input gene IDs are fail to map...

``` r
TCGA_LUAD_gene_map <- TCGA_LUAD_gene_map %>% left_join(TCGA_LUAD_genes)
```

    ## Joining with `by = join_by(ENSEMBL)`

``` r
TCGA_LUSC_gene_map <- TCGA_LUAD_gene_map %>% left_join(TCGA_LUSC_genes)
```

    ## Joining with `by = join_by(ENSEMBL, name, value)`

here we can see some duplicates . By using the map from the cluster
profiler we have now a mapping similar to this :

We can use now this mappings to take our initial expression matrices and
replace the Ensemble ID with the gene Symbols. Important: this will
reduce our set significantly

``` r
old <- dim(TCGA_LUSC_mat)
TCGA_LUSC_mat <- TCGA_LUSC_mat[TCGA_LUSC_gene_map$value,]
new <- dim(TCGA_LUSC_mat)
rownames(TCGA_LUSC_mat) <- TCGA_LUSC_gene_map$SYMBOL

TCGA_LUAD_mat <- TCGA_LUAD_mat[TCGA_LUAD_gene_map$value,]
rownames(TCGA_LUAD_mat) <- TCGA_LUAD_gene_map$SYMBOL
```

from originally to now elements.

## Combination

Now we combine both of them together before doing the comparison

``` r
all.equal(rownames(TCGA_LUAD_mat),rownames(TCGA_LUSC_mat))
```

    ## [1] TRUE

``` r
combined_mat  <- cbind(TCGA_LUAD_mat,TCGA_LUSC_mat)
# The only problem is that we have still not real sample data but the TCGA case numbers

TCGA_lung_meta <- data.frame(cancer_type = c(rep( "LUSC", ncol(TCGA_LUSC_mat)), 
                   rep("LUAD", ncol(TCGA_LUAD_mat))))
```

## PCA analysis

In general for these type of analysis we often dont want to take all but
maybe only the top 2000 most variable ones.

### Raw counts

``` r
most_var_subset_ix <- order(rowVars(combined_mat),decreasing=TRUE)[1:2000]
most_var_subset    <- combined_mat[most_var_subset_ix,]
most_var_pca      <- prcomp(t(most_var_subset),scale.=TRUE)
autoplot(most_var_pca, data=TCGA_lung_meta, color ="cancer_type") + ggtitle("PCA combined lung cancer")
```

![](explore1_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

They are not clearly separable by PC1 which accounts for 20 %
variability, but then they are coming from different experiments and we
took simply the counts. These are though strongly affected by the number
of sequenced reads and not normalized within and neither between
experiments. Now if we assume that this is true, then we would expect
that the sum of counts for a given sample would directly correlate with
it’s PCA value. PC2 is looking better for separation

``` r
nreads <- colSums(most_var_subset)
cor(most_var_pca$x[,1],nreads)
```

    ## [1] 0.9489841

``` r
PCA_overview <- bind_rows(lapply(1:5,function(x){
    data.frame(PC=x,correlation=cor(most_var_pca$x[,x],nreads))
}))
```

So, the first PC is pretty useless because it correlates extremely with
number of counts per sample But the PC2 seems not strongly affected by
the counts of reads.

``` r
TCGA_lung_meta$seq_depth<-nreads
autoplot(most_var_pca, data=TCGA_lung_meta, color ="seq_depth") + scale_color_viridis_b() + ggtitle("PCA combined lung cancer with depth")
```

![](explore1_files/figure-gfm/unnamed-chunk-12-1.png)<!-- --> Yeah,
looks quite obvious indeed.

### CPM normalized

So normally we always normalize, lets do that next.

``` r
combined_mat_cpm <- edgeR::cpm(combined_mat,prior.count=TRUE,log=TRUE)
most_var_cpm_subset_ix <- order(rowVars(combined_mat_cpm),decreasing=TRUE)[1:2000]
most_var_cpm_subset    <- combined_mat_cpm[most_var_cpm_subset_ix,]
most_var_cpm_pca      <- prcomp(t(most_var_cpm_subset),scale.=TRUE)
autoplot(most_var_cpm_pca, data=TCGA_lung_meta, color ="cancer_type") + ggtitle("PCA log2CPM combined lung cancer")
```

![](explore1_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Yeah, that looks much more convincing obviously.

## Heatmap

Here we choose again our 2000 most variable genes

``` r
annot_lung <- HeatmapAnnotation(df = TCGA_lung_meta, col= list(cancer_type= c("LUSC"="red","LUAD"="steelblue")))
Heatmap(most_var_cpm_subset,
    name="Lung cancer log2CPM",
    top_annotation = annot_lung,
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    border = TRUE,
    row_km = 6
)
```

![](explore1_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

So we can see very nicely now the cluster of genes which are up in both,
down in both or really different between both cancer type annotations.
The samples are as well predominantly well separating.

# In depth analysis

Now since we see major differences we can check if some of the know
markers are popping up here as well.

> LUAD originates in the alveolar epithelial cells and is often
> associated with EGFR mutations, particularly in non-smokers

> NKX2-1 (TTF-1) Lung epithelial differentiation marker. Highly
> expressed in LUAD but absent in LUSC (TCGA, IHC studies).

> NAPSA (Napsin A). Aspartic protease in surfactant protein processing.
> High specificity for LUAD (IHC studies).

The first one is a bit more tricky to check as we will have to circle
back to the mutation data. Lets do that later, first have a look at the
expression of these genes and whether they stand out somehow.

``` r
TCGA_lung_meta$NAPSA <- combined_mat_cpm["NAPSA", ]
TCGA_lung_meta$TFF1 <- combined_mat_cpm["TFF1", ]
TCGA_lung_meta$EGFR <- combined_mat_cpm["EGFR", ]

autoplot(most_var_cpm_pca, data = TCGA_lung_meta , color ="TFF1") +
        scale_color_viridis_b() +
        facet_wrap(~ cancer_type) + 
        ggtitle("TCGA TFF1")
```

![](explore1_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
autoplot(most_var_cpm_pca, data = TCGA_lung_meta , color ="NAPSA") +
        scale_color_viridis_b() +
        facet_wrap(~ cancer_type) + 
        ggtitle("TCGA NAPSA")
```

![](explore1_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

``` r
autoplot(most_var_cpm_pca, data = TCGA_lung_meta , color ="EGFR") +
        scale_color_viridis_b() +
        facet_wrap(~ cancer_type) + 
        ggtitle("TCGA EGFR")
```

![](explore1_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->

They are obviously the identical plots, and only the color labels are
changing. Lets have a quick look at the distribution of these:

``` r
ggplot(TCGA_lung_meta,aes(x=TFF1, fill=cancer_type)) + facet_wrap(~cancer_type)+ geom_density() + ggtitle("TFF1 expression")
```

![](explore1_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
ggplot(TCGA_lung_meta,aes(x=NAPSA, fill=cancer_type)) + facet_wrap(~cancer_type)+ geom_density() + ggtitle("NAPSA expression")
```

![](explore1_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

``` r
ggplot(TCGA_lung_meta,aes(x=EGFR, fill=cancer_type)) + facet_wrap(~cancer_type)+ geom_density() + ggtitle("EGFR expression")
```

![](explore1_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->

TTF1 does on both plots look not very impressive - meaning already low
expression and similar distribution. Same for the EGFR expression. But
NAPSA looks indeed already quite different between both types. It seems
that the PC1 is pretty well separating the 2 cancer types based on the
expression of NAPSA, except for some samples. The question might be if
these are mislabeled ?

There are other potential markers that are known in the literature:

> TP63 (ΔNp63). Transcription factor crucial for squamous
> differentiation.High expression in LUSC, absent in LUAD (TCGA, IHC
> studies).

> KRT5 (Cytokeratin 5). Cytoskeletal protein in basal epithelial cells.
> Expressed in LUSC, absent in LUAD.

Lets do the same game as for the others before:

``` r
TCGA_lung_meta$TP63 <- combined_mat_cpm["TP63", ]
TCGA_lung_meta$KRT5 <- combined_mat_cpm["KRT5", ]

autoplot(most_var_cpm_pca, data = TCGA_lung_meta , color ="TP63") +
        scale_color_viridis_b() +
        facet_wrap(~ cancer_type) + 
        ggtitle("TCGA TP63")
```

![](explore1_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
autoplot(most_var_cpm_pca, data = TCGA_lung_meta , color ="KRT5") +
        scale_color_viridis_b() +
        facet_wrap(~ cancer_type) + 
        ggtitle("TCGA KRT5")
```

![](explore1_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

Indeed, these ones are crazy specific.

``` r
ggplot(TCGA_lung_meta,aes(x=TP63, fill=cancer_type)) + facet_wrap(~cancer_type)+ geom_density() + ggtitle("TP63 expression")
```

![](explore1_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
ggplot(TCGA_lung_meta,aes(x=KRT5, fill=cancer_type)) + facet_wrap(~cancer_type)+ geom_density() + ggtitle("KRT5 expression")
```

![](explore1_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

And we see again the same effect that there are some LUAD ones which
seem bizzare and in the PC1 below 0. Lets take these samples and have a
better look at them

``` r
TCGA_lung_meta$PC1_low <- most_var_cpm_pca$x[, 1] < 0 
TCGA_lung_meta$PC1_high <- most_var_cpm_pca$x[, 1] > 0
# the ones which are lower 0 but annotated LUAD
```

How many samples are these which are LUAD annotated but PC1 \< 0 : 115
That’s still many and unlikely to be a labeling problem as it represents
in totla 19.1666667 percent of the samples How many samples are these
which are LUSC annotated but PC1 \> 0 : 45 The other way around we have
only 8.0071174 percent of the samples. That is more likely to be a
mislabeling - but I doubt it. More likely is that it is not a perfect
marker or our PC is not separating them perfectly.

# SNV comparison

As mentioned earlier :

> LUAD originates in the alveolar epithelial cells and is often
> associated with EGFR mutations, particularly in non-smokers

Lets have a look first at this marker

``` r
# This one should only be done once and the rds file been kept
# Therefore it is excluded from the report generation
LUSC_query_snv <- GDCquery(project = "TCGA-LUSC",
                    data.category = "Simple Nucleotide Variation",
                    access = "open",
                    data.type = "Masked Somatic Mutation", 
                    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")

GDCdownload(LUSC_query_snv)
TCGA_LUSC_snv_data <- GDCprepare(LUSC_query_snv)
saveRDS(TCGA_LUSC_snv_data, "./data/TCGA_LUSC_SNV_SummarizedExperiment.rds")
```

``` r
# This one should only be done once and the rds file been kept
# Therefore it is excluded from the report generation
LUAD_query_snv <- GDCquery(project = "TCGA-LUAD",
                    data.category = "Simple Nucleotide Variation",
                    access = "open",
                    data.type = "Masked Somatic Mutation", 
                    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")

GDCdownload(LUAD_query_snv)
TCGA_LUAD_snv_data <- GDCprepare(LUAD_query_snv)
saveRDS(TCGA_LUAD_snv_data, "./data/TCGA_LUAD_SNV_SummarizedExperiment.rds")
```

``` r
maf <- getMC3MAF()
```

``` r
LUAD__snv_data<- readRDS("./data/TCGA_LUAD_SNV_SummarizedExperiment.rds")
LUSC__snv_data<- readRDS("./data/TCGA_LUSC_SNV_SummarizedExperiment.rds")
# save the raw counts matrix 
```

# Survival modelling

For this part, I am switching to this nice
[tutorial](https://ocbe-uio.github.io/survomics/survomics.html) which
has been published alongside this nice review [Tutorial on survival
modeling with applications to omics
data](https://academic.oup.com/bioinformatics/article/40/3/btae132/7623091).

Essentially we are still remaining with the same cancer types but will
now add the survival data. Unfortunately, I am running immmediately here
into some known
[issues](https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/639).
Therefore, currently a need to install devtools and directly from the
git-repository….

``` r
# taken from here https://stackoverflow.com/questions/24519794/r-max-function-ignore-na
my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
cancer_types <- c("TCGA-LUAD","TCGA-LUSC")

clin <- NULL
tmp  <- TCGAbiolinks::GDCquery_clinic(project = "TCGA-LUSC", type = "clinical")
tmp2 <-TCGAbiolinks::GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
clin <- rbind(clin, tmp[, c(
        "project", "submitter_id", "vital_status",
        "days_to_last_follow_up", "days_to_death",
        "age_at_diagnosis", "gender", "race",
        "ethnicity", "ajcc_pathologic_t"
)])

clin <- rbind(clin, tmp2[, c(
        "project", "submitter_id", "vital_status",
        "days_to_last_follow_up", "days_to_death",
        "age_at_diagnosis", "gender", "race",
        "ethnicity", "ajcc_pathologic_t"
)])
```

so now we have for the 2 Lung cancer the following clinical data
imported : project, submitter_id, vital_status, days_to_last_follow_up,
days_to_death, age_at_diagnosis, gender, race, ethnicity,
ajcc_pathologic_t Next we want to get the days to death and days to the
last follow up, lets have a look how many complete, missing and type of
information we have For - time to death: 63.3608815 are `NA` and
missing - time to last follow up: 99.9081726 are `NA` and missing

Next we are converting the resulting times into years instead of days,
same as well for the age at diagnosis. Then we use the maximum of both
values as the endpoint.

``` r
# convert into years
# I got confused because I had many inf values. Turns out that the standard max function returns
# inf if both values are actually "NA" --> which I do not like :)
clin$time <- apply(clin[,c("days_to_death", "days_to_last_follow_up")],1,my.max)/365
clin$age_at_diagnosis<-clin$age_at_diagnosis/365
clin <- clin[, c("project", "submitter_id", "vital_status", "time", "gender", "age_at_diagnosis", "race", "ethnicity")]

ggplot(clin, aes(x=age_at_diagnosis, fill=gender)) + geom_histogram() + ggtitle("Age at diagnosis") + facet_wrap(~project)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 120 rows containing non-finite outside the scale range
    ## (`stat_bin()`).

![](explore1_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
ggplot(clin, aes(x=time, fill=gender)) + geom_histogram() + ggtitle("survival time") + facet_wrap(~project)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 689 rows containing non-finite outside the scale range
    ## (`stat_bin()`).

![](explore1_files/figure-gfm/unnamed-chunk-25-2.png)<!-- -->

``` r
# extract patients with positive overall survival time
clin <- clin[(!is.na(clin$time )) & (clin$vital_status %in% c("Alive", "Dead")), ]
```

Now we only extracted patients with positive overall surivival time and
have in total 400 observations.

``` r
table1 <- clin %>%
  dplyr::count(vital_status, gender,project, ethnicity) %>%
  group_by(vital_status) %>%
  mutate(prop = prop.table(n))
print(table1)
```

    ## # A tibble: 16 × 6
    ## # Groups:   vital_status [2]
    ##    vital_status gender project   ethnicity                  n    prop
    ##    <chr>        <chr>  <chr>     <chr>                  <int>   <dbl>
    ##  1 Alive        male   TCGA-LUSC not reported               1 1      
    ##  2 Dead         female TCGA-LUAD Unknown                    4 0.0100 
    ##  3 Dead         female TCGA-LUAD hispanic or latino         1 0.00251
    ##  4 Dead         female TCGA-LUAD not hispanic or latino    69 0.173  
    ##  5 Dead         female TCGA-LUAD not reported              22 0.0551 
    ##  6 Dead         female TCGA-LUSC hispanic or latino         1 0.00251
    ##  7 Dead         female TCGA-LUSC not hispanic or latino    29 0.0727 
    ##  8 Dead         female TCGA-LUSC not reported              20 0.0501 
    ##  9 Dead         male   TCGA-LUAD Unknown                    4 0.0100 
    ## 10 Dead         male   TCGA-LUAD hispanic or latino         2 0.00501
    ## 11 Dead         male   TCGA-LUAD not hispanic or latino    58 0.145  
    ## 12 Dead         male   TCGA-LUAD not reported              24 0.0602 
    ## 13 Dead         male   TCGA-LUSC Unknown                    1 0.00251
    ## 14 Dead         male   TCGA-LUSC hispanic or latino         5 0.0125 
    ## 15 Dead         male   TCGA-LUSC not hispanic or latino    96 0.241  
    ## 16 Dead         male   TCGA-LUSC not reported              63 0.158

``` r
ggplot(clin, aes(x=time, fill=vital_status)) + geom_histogram() + ggtitle("survival time") + facet_wrap(~project)
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](explore1_files/figure-gfm/unnamed-chunk-26-1.png)<!-- --> I think
the histogram is much better suited for the visualization than the
lollipop they show ![image in their
repo](https://ocbe-uio.github.io/survomics/fig/TCGA_survival.png)

Now we take our perviously normalized RNAseq data and do survival
analysis using as well the clinical/demographic variables.

### LUSC

I added to the typical values as well “\$tobacco_smoking_status” to
inform us about the most likely source of lung cancer.

``` r
# we are reusing the followng elements again
#TCGA_lung_meta
#combined_mat_cpm
#rownames(TCGA_lung_meta) <- colnames(combined_mat_cpm)

meta_LUSC <- colData(TCGA_LUSC_data)[,c("project_id", "submitter_id", "age_at_diagnosis", 
        "ethnicity", "gender", "days_to_death", 
        "days_to_last_follow_up", "vital_status","treatments","tobacco_smoking_status","paper_T.stage")]

# now we only differentiate if/or not treatment
meta_LUSC$treatments <- unlist(lapply(meta_LUSC$treatments, function(y) {
  any(y$treatment_or_therapy == "yes")
}))

smoker = c("Current Reformed Smoker for < or = 15 yrs","Current Smoker","Current Reformed Smoker, Duration Not Specified")
no_smoker = c("Lifelong Non-Smoker","Current Reformed Smoker for > 15 yrs")

# same for smoker
meta_LUSC$tobacco_smoking_status[meta_LUSC$tobacco_smoking_status=="Unknown" | meta_LUSC$tobacco_smoking_status=="Not Reported"] <- NA
meta_LUSC$tobacco_smoking_status[meta_LUSC$tobacco_smoking_status %in% smoker ] <- TRUE
meta_LUSC$tobacco_smoking_status[meta_LUSC$tobacco_smoking_status %in% no_smoker] <- FALSE

# same for tumor stage, it is already very sparse and then
# we have sub-groups. Lets try to make just 4 groups instead
t1 <- c("T1","T1a","T1b")
t2 <- c("T2","T2a","T2b")
meta_LUSC$paper_T.stage[meta_LUSC$paper_T.stage %in% t1] <- "T1"
meta_LUSC$paper_T.stage[meta_LUSC$paper_T.stage %in% t2] <- "T2"
meta_LUSC$paper_T.stage <- droplevels(meta_LUSC$paper_T.stage)
# now I have already edge and dont want to install DESEQ2 just for normalization
# therefore we do simply here a cpm transformation
# thanks to the above steps, we already have gene names instead of identifiers

LUSC_RNA_cpm <- edgeR::cpm(TCGA_LUSC_mat,prior.count=TRUE,log=TRUE)
# same here as above, using our own function for double NAs
meta_LUSC$time <- apply(meta_LUSC[, c("days_to_death", "days_to_last_follow_up")], 1, my.max) / 365
meta_LUSC$status <- meta_LUSC$vital_status
meta_LUSC$age <- meta_LUSC$age_at_diagnosis / 365
clin_LUSC <- meta_LUSC
#clin_LUSC <- subset(meta_LUSC, !duplicated(submitter_id) & !is.na(time) )
clin_LUSC <- clin_LUSC[order(clin_LUSC$submitter_id), ]
clin_LUSC <- as.data.frame(clin_LUSC)

LUSC_RNA_cpm <- LUSC_RNA_cpm[, rownames(clin_LUSC)]
```

By removing duplicated submitter ID, empty estimation times and unknown
ages, we reduce the set of patients from originally 562 to 562.

#### Plotting expression

So normally we always normalize, lets do that next.

``` r
most_var_cpm_subset_ix <- order(rowVars(LUSC_RNA_cpm),decreasing=TRUE)[1:5000]
most_var_cpm_subset    <- LUSC_RNA_cpm[most_var_cpm_subset_ix,]
most_var_cpm_pca      <- prcomp(t(most_var_cpm_subset),scale.=TRUE)
autoplot(most_var_cpm_pca, data=clin_LUSC, color ="treatments") + ggtitle("PCA log2CPM LUSC subset")
```

![](explore1_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
autoplot(most_var_cpm_pca, data=clin_LUSC, color ="gender") + ggtitle("PCA log2CPM LUSC subset")
```

![](explore1_files/figure-gfm/unnamed-chunk-28-2.png)<!-- -->

``` r
autoplot(most_var_cpm_pca, data=clin_LUSC, color ="ethnicity") + ggtitle("PCA log2CPM LUSC subset")
```

![](explore1_files/figure-gfm/unnamed-chunk-28-3.png)<!-- -->

``` r
autoplot(most_var_cpm_pca, data=clin_LUSC, color ="tobacco_smoking_status") + ggtitle("PCA log2CPM LUSC subset")
```

![](explore1_files/figure-gfm/unnamed-chunk-28-4.png)<!-- -->

``` r
autoplot(most_var_cpm_pca, data=clin_LUSC, color ="paper_T.stage") + ggtitle("PCA log2CPM LUSC subset")
```

![](explore1_files/figure-gfm/unnamed-chunk-28-5.png)<!-- -->

There is nothing that clearly separates the 5000 most variable genes in
the PCA

#### Heatmap

Here we choose again our 2000 most variable genes

``` r
annot_lung <- HeatmapAnnotation(
        gender = clin_LUSC$gender,
        ethnicity=clin_LUSC$ethnicity,
        vital_status=clin_LUSC$vital_status,
        treatmens=clin_LUSC$treatments,
        smoking=clin_LUSC$tobacco_smoking_status,
        tumor_stage=clin_LUSC$paper_T.stage
        )

Heatmap(most_var_cpm_subset,
    name="Lung cancer log2CPM",
    top_annotation = annot_lung,
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    border = TRUE,
    row_km = 6
)
```

    ## `use_raster` is automatically set to TRUE for a matrix with more than
    ## 2000 rows. You can control `use_raster` argument by explicitly setting
    ## TRUE/FALSE to it.
    ## 
    ## Set `ht_opt$message = FALSE` to turn off this message.

    ## 'magick' package is suggested to install to give better rasterization.
    ## 
    ## Set `ht_opt$message = FALSE` to turn off this message.

![](explore1_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

We can see similarly in the heatmap that none of our annotated values is
really helping here in the clustering of genes.

#### Nonparametric surival analysis

Lets plot some Kaplan-Maier survival probability curves

``` r
clin_LUSC$bstatus[clin_LUSC$status == "Dead"] <- 1 
clin_LUSC$bstatus[clin_LUSC$status == "Alive"] <- 0
clin_LUSC$bstatus <- as.numeric(clin_LUSC$bstatus)
LUSC_sfit <- survival::survfit(Surv(time, bstatus) ~ 1, data = clin_LUSC)
LUSC_sfit
```

    ## Call: survfit(formula = Surv(time, bstatus) ~ 1, data = clin_LUSC)
    ## 
    ##    310 observations deleted due to missingness 
    ##        n events median 0.95LCL 0.95UCL
    ## [1,] 252    251   1.49    1.23    1.83

Now we have a fit, maybe important to mention that we have only 306
observations where the patient survived at all Lets now plot a curve

``` r
LUSC_ggsurv <- survminer::ggsurvplot(LUSC_sfit,
  conf.int = TRUE, risk.table = TRUE,
  xlab = "Time since diagnosis (year)",
  legend = "none", surv.median.line = "hv"
)

LUSC_ggsurv$plot <- LUSC_ggsurv$plot + annotate("text", x = 20, y = 0.9, label = "+  Censor")
LUSC_ggsurv
```

![](explore1_files/figure-gfm/km_curve1-1.png)<!-- -->

Not a very great outcome, lets compare treatment vs non-treated.

``` r
LUSC_comp <- survival::survdiff(Surv(time, bstatus) ~ treatments, data = clin_LUSC)
LUSC_sfit2 <- survival::survfit(Surv(time, bstatus) ~ treatments, data = clin_LUSC)
LUSC_ggsurv <- survminer::ggsurvplot(LUSC_sfit2,
  conf.int = TRUE, risk.table = TRUE,
  xlab = "Time since diagnosis (year)", legend = c(.6, .9),
  legend.labs = c("No", "Yes"), legend.title = "Treatment",
  risk.table.y.text.col = TRUE, risk.table.y.text = FALSE, pval=TRUE
)
LUSC_ggsurv
```

![](explore1_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

Hmm not looking particularly promising concerning treatment for this
data-set. The p-value is non-significant, so treatment alone does not
significantly sway prognosis in this dataset. I tried the same with sex,
smoke status and tumor stage, but no luck.

Lets have a look with a Cox model if we can identify any factor which
might actually be interesting here.

``` r
coxph(Surv(time, bstatus) ~ paper_T.stage + treatments, data = clin_LUSC)
```

    ## Call:
    ## coxph(formula = Surv(time, bstatus) ~ paper_T.stage + treatments, 
    ##     data = clin_LUSC)
    ## 
    ##                    coef exp(coef) se(coef)      z       p
    ## paper_T.stageT2  0.6118    1.8438   0.2711  2.256 0.02404
    ## paper_T.stageT3  1.2007    3.3224   0.4528  2.652 0.00801
    ## paper_T.stageT4  0.9998    2.7178   0.5118  1.954 0.05076
    ## treatmentsTRUE  -0.6538    0.5201   0.2432 -2.688 0.00718
    ## 
    ## Likelihood ratio test=14.07  on 4 df, p=0.007065
    ## n= 90, number of events= 90 
    ##    (472 observations deleted due to missingness)

We can see here that the tumor stage at which it was detected has strong
implications on the outcome. Especially if we put them here in the model
in combination with the treatment. Lets plot again survival curves, but
now including as well the stage of the tumor

``` r
LUSC_sfit2 <- survival::survfit(Surv(time, bstatus) ~ treatments + paper_T.stage, data = clin_LUSC)
LUSC_ggsurv <- survminer::ggsurvplot(LUSC_sfit2,
  risk.table = TRUE,
  xlab = "Time since diagnosis (year)",
  risk.table.y.text.col = TRUE, risk.table.y.text = FALSE, pval=TRUE
)
LUSC_ggsurv
```

![](explore1_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

> The Cox model assumes proportional hazards and log-linearity of the
> covariates. To check the log-linearity for a clinical or demographic
> variable, e.g. age, we can fit a penalized smoothing spline for age
> effect

Lets see what that gives . Note: unfortunately, did the provided
functions not tolerate NA values. I need to re-visit this later again.

### LUAD survival analysis

I added to the typical values as well “\$tobacco_smoking_status” to
inform us about the most likely source of lung cancer.

``` r
# we are reusing the followng elements again
#TCGA_lung_meta
#combined_mat_cpm
#rownames(TCGA_lung_meta) <- colnames(combined_mat_cpm)

meta_LUAD <- colData(TCGA_LUAD_data)[,c("project_id", "submitter_id", "age_at_diagnosis", 
        "ethnicity", "gender", "days_to_death", 
         "vital_status","treatments","tobacco_smoking_status","paper_T.stage","paper_N.stage")]



# now we only differentiate if/or not treatment
meta_LUAD$treatments <- unlist(lapply(meta_LUAD$treatments, function(y) {
  any(y$treatment_or_therapy == "yes")
}))

smoker = c("Current Reformed Smoker for < or = 15 yrs","Current Smoker","Current Reformed Smoker, Duration Not Specified")
no_smoker = c("Lifelong Non-Smoker","Current Reformed Smoker for > 15 yrs")

# same for smoker
meta_LUAD$tobacco_smoking_status[meta_LUAD$tobacco_smoking_status=="Unknown" | meta_LUAD$tobacco_smoking_status=="Not Reported"] <- NA
meta_LUAD$tobacco_smoking_status[meta_LUAD$tobacco_smoking_status %in% smoker ] <- TRUE
meta_LUAD$tobacco_smoking_status[meta_LUAD$tobacco_smoking_status %in% no_smoker] <- FALSE

# same for tumor stage, it is already very sparse and then
# we have sub-groups. Lets try to make just 4 groups instead
t1 <- c("T1","T1a","T1b")
t2 <- c("T2","T2a","T2b")
meta_LUAD$paper_T.stage[meta_LUAD$paper_T.stage %in% t1] <- "T1"
meta_LUAD$paper_T.stage[meta_LUAD$paper_T.stage %in% t2] <- "T2"
meta_LUAD$paper_T.stage <- droplevels(meta_LUAD$paper_T.stage)
# now I have already edge and dont want to install DESEQ2 just for normalization
# therefore we do simply here a cpm transformation
# thanks to the above steps, we already have gene names instead of identifiers

LUAD_RNA_cpm <- edgeR::cpm(TCGA_LUAD_mat,prior.count=TRUE,log=TRUE)
# same here as above, using our own function for double NAs
meta_LUAD$time <-meta_LUAD$days_to_death/ 365

meta_LUAD$status <- meta_LUAD$vital_status
meta_LUAD$age <- meta_LUAD$age_at_diagnosis / 365
clin_LUAD <- meta_LUAD
#clin_LUSC <- subset(meta_LUSC, !duplicated(submitter_id) & !is.na(time) )
clin_LUAD <- clin_LUAD[order(clin_LUAD$submitter_id), ]
clin_LUAD <- as.data.frame(clin_LUAD)

LUAD_RNA_cpm <- LUAD_RNA_cpm[, rownames(clin_LUAD)]
```

By removing duplicated submitter ID, empty estimation times and unknown
ages, we reduce the set of patients from originally 600 to 600.

#### Plotting expression

So normally we always normalize, lets do that next.

``` r
most_var_cpm_subset_ix <- order(rowVars(LUAD_RNA_cpm),decreasing=TRUE)[1:5000]
most_var_cpm_subset    <- LUAD_RNA_cpm[most_var_cpm_subset_ix,]
most_var_cpm_pca      <- prcomp(t(most_var_cpm_subset),scale.=TRUE)
autoplot(most_var_cpm_pca, data=clin_LUAD, color ="treatments") + ggtitle("PCA log2CPM LUAD subset")
```

![](explore1_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

``` r
autoplot(most_var_cpm_pca, data=clin_LUAD, color ="gender") + ggtitle("PCA log2CPM LUAD subset")
```

![](explore1_files/figure-gfm/unnamed-chunk-35-2.png)<!-- -->

``` r
autoplot(most_var_cpm_pca, data=clin_LUAD, color ="ethnicity") + ggtitle("PCA log2CPM LUAD subset")
```

![](explore1_files/figure-gfm/unnamed-chunk-35-3.png)<!-- -->

``` r
autoplot(most_var_cpm_pca, data=clin_LUAD, color ="tobacco_smoking_status") + ggtitle("PCA log2CPM LUAD subset")
```

![](explore1_files/figure-gfm/unnamed-chunk-35-4.png)<!-- -->

``` r
autoplot(most_var_cpm_pca, data=clin_LUAD, color ="paper_T.stage") + ggtitle("PCA log2CPM LUAD subset")
```

![](explore1_files/figure-gfm/unnamed-chunk-35-5.png)<!-- -->

There is nothing that clearly separates the 5000 most variable genes in
the PCA

#### Heatmap

Here we choose again our 2000 most variable genes

``` r
annot_lung <- HeatmapAnnotation(
        gender = clin_LUAD$gender,
        ethnicity=clin_LUAD$ethnicity,
        vital_status=clin_LUAD$vital_status,
        treatmens=clin_LUAD$treatments,
        smoking=clin_LUAD$tobacco_smoking_status,
        tumor_stage=clin_LUAD$paper_T.stage,
        M_stage=clin_LUAD$paper_M.stage,
        N_stage=clin_LUAD$paper_N.stage)

Heatmap(most_var_cpm_subset,
    name="Lung cancer log2CPM",
    top_annotation = annot_lung,
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    border = TRUE,
    row_km = 6
)
```

    ## `use_raster` is automatically set to TRUE for a matrix with more than
    ## 2000 rows. You can control `use_raster` argument by explicitly setting
    ## TRUE/FALSE to it.
    ## 
    ## Set `ht_opt$message = FALSE` to turn off this message.

    ## 'magick' package is suggested to install to give better rasterization.
    ## 
    ## Set `ht_opt$message = FALSE` to turn off this message.

![](explore1_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

We can see similarly in the heatmap that none of our annotated values is
really helping here in the clustering of genes, neither.

#### Nonparametric surival analysis

Lets plot some Kaplan-Maier survival probability curves

``` r
clin_LUAD$bstatus[clin_LUAD$status == "Dead"] <- 1 
clin_LUAD$bstatus[clin_LUAD$status == "Alive"] <- 0
clin_LUAD$bstatus <- as.numeric(clin_LUAD$bstatus)
LUAD_sfit <- survival::survfit(Surv(time, bstatus) ~ 1, data = clin_LUAD)
LUAD_sfit
```

    ## Call: survfit(formula = Surv(time, bstatus) ~ 1, data = clin_LUAD)
    ## 
    ##    385 observations deleted due to missingness 
    ##        n events median 0.95LCL 0.95UCL
    ## [1,] 215    215   1.71    1.37    2.02

Now we have a fit, maybe important to mention that we have only 381
observations where the patient survived at all Lets now plot a curve

``` r
LUAD_ggsurv <- survminer::ggsurvplot(LUAD_sfit,
  conf.int = TRUE, risk.table = TRUE,
  xlab = "Time since diagnosis (year)",
  legend = "none", surv.median.line = "hv"
)

LUAD_ggsurv$plot <- LUAD_ggsurv$plot + annotate("text", x = 20, y = 0.9, label = "+  Censor")
LUAD_ggsurv
```

![](explore1_files/figure-gfm/km_curve-1.png)<!-- -->

Not a very great outcome, lets compare treatment vs non-treated.

``` r
LUAD_comp <- survival::survdiff(Surv(time, bstatus) ~ treatments, data = clin_LUAD)
LUAD_sfit2 <- survival::survfit(Surv(time, bstatus) ~ treatments, data = clin_LUAD)
LUAD_ggsurv <- survminer::ggsurvplot(LUAD_sfit2,
  conf.int = TRUE, risk.table = TRUE,
  xlab = "Time since diagnosis (year)", legend = c(.6, .9),
  legend.labs = c("No", "Yes"), legend.title = "Treatment",
  risk.table.y.text.col = TRUE, risk.table.y.text = FALSE, pval=TRUE
)

LUAD_ggsurv
```

![](explore1_files/figure-gfm/surv_plot_LUAD1-1.png)<!-- -->

Same game now as above for the LUSC, lets see if we can find factors
which would help to discriminate something.

``` r
coxph(Surv(time, bstatus) ~ paper_N.stage, data = clin_LUAD)
```

    ## Call:
    ## coxph(formula = Surv(time, bstatus) ~ paper_N.stage, data = clin_LUAD)
    ## 
    ##                                coef exp(coef) se(coef)     z       p
    ## paper_N.stageN1              0.3991    1.4905   0.2506 1.592 0.11130
    ## paper_N.stageN2              0.7174    2.0491   0.2709 2.648 0.00809
    ## paper_N.stageN3                  NA        NA   0.0000    NA      NA
    ## paper_N.stage[Not Available]     NA        NA   0.0000    NA      NA
    ## paper_N.stageNX              0.7582    2.1344   0.7328 1.035 0.30088
    ## 
    ## Likelihood ratio test=7.45  on 3 df, p=0.05882
    ## n= 95, number of events= 95 
    ##    (505 observations deleted due to missingness)

I found one element, which seems interesting which is the N-stage from
the paper. This is the extend of nearby affected lymph nodes.

``` r
LUAD_sfit2 <- survival::survfit(Surv(time, bstatus) ~ paper_N.stage, data = clin_LUAD)
LUAD_ggsurv <- survminer::ggsurvplot(LUAD_sfit2,
  risk.table = TRUE,
  xlab = "Time since diagnosis (year)",
  risk.table.y.text.col = TRUE, risk.table.y.text = FALSE, pval=TRUE
)
LUAD_ggsurv
```

![](explore1_files/figure-gfm/surv_plot_LUAD2-1.png)<!-- -->

Here we can see quite nicely how the different grades of lymph node
invasion influence the outcome.
