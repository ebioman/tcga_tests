FROM r-base:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends libxml2-dev
# basic packages
RUN Rscript -e 'install.packages(c("languageserver","BiocManager","dplyr","here"))'
RUN apt-get update && apt-get install -y --no-install-recommends curl libcurl4-openssl-dev
# install BioC and TCGA packages
RUN Rscript -e 'BiocManager::install(c("TCGAbiolinks","SummarizedExperiment"))'
# add if needed small additional packages
RUN Rscript -e 'install.packages(c("janitor","rmarkdown","scatterpie","ggfortify"))'
RUN apt-get update && apt-get install -y --no-install-recommends libfontconfig1-dev
RUN Rscript -e 'BiocManager::install(c("org.Hs.eg.db","ComplexHeatmap","clusterProfiler","edgeR"))'
RUN apt-get update && apt-get install -y --no-install-recommends pandoc
RUN Rscript -e 'install.packages(c("languageserver","rmarkdown","httpgd"))'

