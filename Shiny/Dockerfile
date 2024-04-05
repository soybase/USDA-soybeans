# shiny-verse instead of shiny for devtools & dplyr
FROM rstudio/r-base:3.5.3-focal AS base

RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/focal/2019-04-25"))' > /usr/lib/R/etc/Rprofile.site

# create directory for logs & bookmarks
WORKDIR /srv
RUN rm -rf *

COPY CNV/ CNV/
COPY Genealogy/ Genealogy/
COPY SNP/ SNP/

FROM base AS docs

RUN Rscript -e 'install.packages(c("knitrBootstrap", "xtable"))'

COPY ./statistics.csl ./Refs.bib ./

RUN for dir in CNV Genealogy SNP; do \
  ( cd ${dir} \
  && ln -s ../statistics.csl ../Refs.bib . \
  && Rscript -e 'library(knitrBootstrap); library(rmarkdown); render("Documentation.Rmd", "knitrBootstrap::bootstrap_document")' \
  && rm statistics.csl Refs.bib); done

FROM docker.io/rhub/r-minimal:3.6.3 AS final

# MASS@7.3-58.3: ggplot2 dependency; MASS > 7.3-58.3 needs R >= 4.0
# lattice@0.20-41: Matrix dependency; lattice > 0.20-41 needs R >= 4.0
# lattice@0.21-7: ggplot2 dependency
RUN installr -a 'cairo font-liberation libcurl' -t 'cairo-dev curl-dev gfortran zlib-dev' -d MASS@7.3-58.3 lattice@0.20-41 Cairo ggenealogy ggplot2 knitr magrittr shiny dicook/phyViz@18595fb83c8f604e896408fbe201f1c7f700ab20

WORKDIR /srv

COPY CNV/ CNV/
COPY Genealogy/ Genealogy/
COPY SNP/ SNP/

COPY --from=docs /srv/CNV/Documentation.html CNV/
COPY --from=docs /srv/Genealogy/Documentation.html Genealogy/
COPY --from=docs /srv/SNP/Documentation.html SNP/

RUN find . -name '*.tar.xz' | while read -r f; do tar -C $(dirname "$f") -xJf "$f" && rm "$f"; done

USER nobody
