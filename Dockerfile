FROM 579831337053.dkr.ecr.us-west-2.amazonaws.com/irx-r-base:latest

ENV _R_CHECK_TESTS_NLINES_=0

WORKDIR /src
RUN apt-get -y update
RUN apt-get -y install git curl
RUN git clone https://github.com/metrumresearchgroup/Torsten.git
WORKDIR /src/Torsten/cmdstan
RUN make build

COPY ./ /src/PKPDposterior
WORKDIR /src/PKPDposterior
RUN Rscript -e "devtools::install()"
RUN Rscript -e "install.packages(c('mockery'), repos = '${RSPM_SNAPSHOT}')"
ENV CMDSTAN="/src/Torsten/cmdstan"
