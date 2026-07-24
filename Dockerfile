FROM bioconductor/bioconductor_docker:3.21

RUN Rscript -e "\
  BiocManager::install(c( \
    'VariantAnnotation', \
    'BSgenome', \
    'BSgenome.Hsapiens.UCSC.hg38', \
    'GenomeInfoDb', \
    'Rsamtools', \
    'GenomicRanges', \
    'IRanges', \
    'S4Vectors', \
    'SummarizedExperiment', \
    'Biostrings', \
    'DECIPHER', \
    'biomaRt' \
  ), ask=FALSE, update=FALSE)"

RUN Rscript -e "\
  BiocManager::install('maftools', ask=FALSE, update=FALSE)"

RUN Rscript -e "\
  install.packages(c('UniprotR', 'httr', 'stringr', 'rlang', \
    'knitr', 'kableExtra', 'ggplot2', 'dplyr', 'quarto'), \
    repos='https://cloud.r-project.org')"

RUN Rscript -e "\
  remotes::install_git( \
    'https://gitlab.uvm.edu/vigr/projects/vigr/BSRBioinformaticsTools.git', \
    upgrade='never')"

COPY analysis.R /opt/analysis.R
COPY helperFunctions.R /opt/helperFunctions.R
COPY template.qmd /opt/template.qmd

ENTRYPOINT ["Rscript", "/opt/analysis.R"]
