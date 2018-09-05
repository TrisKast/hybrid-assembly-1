FROM continuumio/miniconda

COPY environment_reduced.yml /
RUN conda env create -f /environment_reduced.yml && conda clean -a
ENV PATH /opt/conda/envs/assembly-env/bin:$PATH

# Install MaSuRCA 3.2.7
RUN apt-get update && apt-get install -y g++ libboost-all-dev zlib1g-dev libbz2-dev make
RUN curl -fsSL https://github.com/alekseyzimin/masurca/raw/master/MaSuRCA-3.2.8.tar.gz -o /opt/MaSuRCA-3.2.8.tar.gz
RUN cd /opt/; tar -xzvf MaSuRCA-3.2.8.tar.gz; cd MaSuRCA-3.2.8; ./install.sh
ENV PATH $PATH:/opt/MaSuRCA-3.2.8/bin

RUN conda install -c bioconda samtools

#Environment for nanoqc
COPY nanoqc-env.yml /
RUN conda env create -f /nanoqc-env.yml

RUN git clone https://github.com/MariaNattestad/Assemblytics.git
ENV PATH $PATH:/Assemblytics

RUN conda create --name python_env python=2.7

#Environment for AsmVar
#COPY asmvar-env.yml /
#RUN conda env create -f /asmvar-env.yml

#Install AsmVar and activate base environment afterwards again
#RUN git clone https://github.com/bioinformatics-centre/AsmVar.git
#RUN cd AsmVar/src/AsmvarDetect; make




# minikraken DB 
# RUN mkdir /kraken_db/ && cd /kraken_db/ && wget https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz && tar xf minikraken.tgz && rm minikraken.tgz
# ENV KRAKEN_DB_PATH="/kraken_db:${KRAKEN_DB_PATH}"
