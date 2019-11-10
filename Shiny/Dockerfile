# shiny-verse instead of shiny for devtools & dplyr
FROM rocker/shiny-verse:3.5.3

RUN installGithub.r --update=FALSE dicook/phyViz
# https://github.com/rocker-org/shiny/issues/68
RUN install2.r --error knitrBootstrap \
  && rm /tmp/downloaded_packages/*

# create directory for logs & bookmarks
WORKDIR /srv/shiny-server
RUN rm -rf *

COPY ./statistics.csl ./Refs.bib ./
COPY Genealogy/ Genealogy/
RUN cd Genealogy \
  && Rscript -e 'library(knitrBootstrap); library(rmarkdown); render("Documentation.Rmd", "knitrBootstrap::bootstrap_document")'

COPY CNV/ CNV/
RUN cd CNV \
  && Rscript -e 'library(knitrBootstrap); library(rmarkdown); render("Documentation.Rmd", "knitrBootstrap::bootstrap_document")'

RUN rm -f ./statistics.csl ./Refs.bib

EXPOSE 3838
