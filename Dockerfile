FROM docker.pkg.github.com/califano-lab/phosphoviper.db/phosphoviper.db:latest

# install phosphoviper
ADD ./ ./
RUN R CMD INSTALL ./

# docker build -t docker.pkg.github.com/califano-lab/phosphoviper/phosphoviper:latest ./
# docker push docker.pkg.github.com/califano-lab/phosphoviper/phosphoviper:latest
