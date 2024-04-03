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

FROM base AS final

RUN Rscript -e 'install.packages(c("animint", "devtools", "ggenealogy", "ggplot2", "knitr", "magrittr", "shiny"))'
RUN Rscript -e 'library(devtools); install_github("dicook/phyViz@18595fb83c8f604e896408fbe201f1c7f700ab20")'

COPY --from=docs /srv/CNV/Documentation.html CNV/
COPY --from=docs /srv/Genealogy/Documentation.html Genealogy/
COPY --from=docs /srv/SNP/Documentation.html SNP/

# libglpk40 needed by igraph
RUN apt update \
  && apt install -y --no-install-recommends libglpk40 xz-utils \
  && find . -name '*.tar.xz' | while read -r f; do tar -C $(dirname "$f") -xJf "$f" && rm "$f"; done \
  && apt remove -y xz-utils \
  && apt autoremove \
  && rm -rf /var/lib/apt/lists/*