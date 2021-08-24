FROM docker.pkg.github.com/califano-lab/vespa.db/vespa.db:latest

# install vespa
ADD ./ ./
RUN R CMD INSTALL ./

# docker build -t docker.pkg.github.com/califano-lab/vespa/vespa:latest ./
# docker push docker.pkg.github.com/califano-lab/vespa/vespa:latest
