# shiny-verse instead of shiny for devtools & dplyr
FROM rocker/shiny-verse:3.5.3

RUN installGithub.r --update=FALSE dicook/phyViz@18595fb83c8f604e896408fbe201f1c7f700ab20
# https://github.com/rocker-org/shiny/issues/68
RUN install2.r --error --skipinstalled \
  ggenealogy \
  knitrBootstrap \
  && rm /tmp/downloaded_packages/*

# create directory for logs & bookmarks
WORKDIR /srv/shiny-server
RUN rm -rf *

COPY CNV/ CNV/
COPY Genealogy/ Genealogy/
COPY SNP/ SNP/

RUN find . -name '*.tar.xz' | while read -r f; do tar -C $(dirname "$f") -xJf "$f"; done

COPY ./statistics.csl ./Refs.bib ./
RUN for dir in CNV Genealogy SNP; do \
  ( cd ${dir} \
  && ln -s ../statistics.csl ../Refs.bib . \
  && Rscript -e 'library(knitrBootstrap); library(rmarkdown); render("Documentation.Rmd", "knitrBootstrap::bootstrap_document")' \
  && rm statistics.csl Refs.bib); done

RUN rm -f ./statistics.csl ./Refs.bib
