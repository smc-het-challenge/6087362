FROM r-base

# install packages
RUN echo 'install.packages("matrixStats", repos = "http://cran.us.r-project.org", dependencies = TRUE)' > /tmp/packages.R && Rscript /tmp/packages.R
