FROM continuumio/miniconda

COPY base-env.yml /
RUN conda env create -f /base-env.yml && conda clean -a

# Install MaSuRCA 3.3.1
RUN apt-get update && apt-get install -y g++ libboost-all-dev zlib1g-dev libbz2-dev make
RUN curl -fsSL https://github.com/alekseyzimin/masurca/raw/master/MaSuRCA-3.3.1.tar.gz -o /opt/MaSuRCA-3.3.1.tar.gz
RUN cd /opt/; tar -xzvf MaSuRCA-3.3.1.tar.gz; cd MaSuRCA-3.3.1; ./install.sh
ENV PATH $PATH:/opt/MaSuRCA-3.3.1/bin

#Install samtools manually because of dependencies issues
RUN conda install -c bioconda -c r samtools --override-channels

#Environment for nanoqc
COPY nanoqc-env.yml /
RUN conda env create -f /nanoqc-env.yml

#Install Assemblytics
RUN git clone https://github.com/MariaNattestad/Assemblytics.git
ENV PATH $PATH:/Assemblytics

#Create python2.7 environment for Assemblytics
COPY python2_7-env.yml /
RUN conda env create -f /python2_7-env.yml

ENV PATH /opt/conda/envs/assembly-env/bin:$PATH


#minikraken DB
#RUN mkdir /kraken_db/ && cd /kraken_db/ && wget https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz && tar xf minikraken.tgz && rm minikraken.tgz
#ENV KRAKEN_DB_PATH="/kraken_db:${KRAKEN_DB_PATH}"
