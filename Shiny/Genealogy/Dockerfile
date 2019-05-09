# shiny-verse instead of shiny for devtools & dplyr
FROM rocker/shiny-verse:3.6.1

RUN installGithub.r --update=FALSE dicook/phyViz
# https://github.com/rocker-org/shiny/issues/68
RUN install2.r --error knitrBootstrap \
    && rm /tmp/downloaded_packages/*

# create directory for logs & bookmarks
WORKDIR /srv/shiny-server
COPY Documentation.Rmd SelectGeneology.R server.R ui.R soybeanGenealogy.csv.gz ./
RUN Rscript -e 'library(knitrBootstrap); library(rmarkdown); render("Documentation.Rmd", "knitrBootstrap::bootstrap_document")'

EXPOSE 3838
