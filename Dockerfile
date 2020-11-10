FROM r-base

# install Debian dependencies
RUN apt-get update && apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev

# install R dependencies
RUN R -e "install.packages(c('pROC','data.table','pbapply','seqinr','stringr','mixtools','plyr','reshape','tidyr','BiocManager') ,dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('viper','BatchQC','preprocessCore','sva','limma'))"

# install phosphoVIPER
ADD ./ ./
RUN R CMD INSTALL ./

# docker build -t docker.pkg.github.com/califano-lab/phosphoviper/phosphoviper:latest ./
# docker push docker.pkg.github.com/califano-lab/phosphoviper/phosphoviper:latest
