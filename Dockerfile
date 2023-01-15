FROM ghcr.io/califano-lab/vespa.db/vespa.db:latest

# install vespa
ADD ./ ./
RUN R CMD INSTALL ./

# docker build -t ghcr.io/califano-lab/vespa/vespa:latest ./
# docker push ghcr.io/califano-lab/vespa/vespa:latest
