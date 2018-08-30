FROM continuumio/miniconda

COPY environment_reduced.yml /
RUN conda env create -f /environment_reduced.yml && conda clean -a
ENV PATH /opt/conda/envs/assembly-env/bin:$PATH

# Install MaSuRCA 3.2.7
RUN apt-get update && apt-get install -y g++ libboost-all-dev zlib1g-dev libbz2-dev make
RUN curl -fsSL https://github.com/alekseyzimin/masurca/raw/master/MaS-c bioconda nanoqcuRCA-3.2.7.tar.gz -o /opt/MaSuRCA-3.2.7.tar.gz
RUN cd /opt/; tar -xzvf MaSuRCA-3.2.7.tar.gz; cd MaSuRCA-3.2.7; ./install.sh
ENV PATH $PATH:/opt/MaSuRCA-3.2.7/bin

RUN conda install -c bioconda samtools

COPY nanoqc-env.yml /
RUN conda env create -f /nanoqc-env.yml

# minikraken DB 
# RUN mkdir /kraken_db/ && cd /kraken_db/ && wget https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz && tar xf minikraken.tgz && rm minikraken.tgz
# ENV KRAKEN_DB_PATH="/kraken_db:${KRAKEN_DB_PATH}"
