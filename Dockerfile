FROM continuumio/miniconda

COPY environment_reduced.yml /
RUN conda env create -f /environment_reduced.yml && conda clean -a
ENV PATH /opt/conda/envs/assembly-env/bin:$PATH

COPY nanoqc-env.yml /
RUN conda env create -f /nanoqc-env.yml

# Install MaSuRCA 3.2.7
RUN apt-get update && apt-get install -y g++ libboost-all-dev zlib1g-dev libbz2-dev make
RUN curl -fsSL https://github.com/alekseyzimin/masurca/raw/master/MaSuRCA-3.2.7.tar.gz -o /opt/MaSuRCA-3.2.7.tar.gz
RUN cd /opt/; tar -xzvf MaSuRCA-3.2.7.tar.gz; cd MaSuRCA-3.2.7; ./install.sh
ENV PATH $PATH:/opt/MaSuRCA-3.2.7/bin

#Install newest version of QUASt for the "large genome" option
RUN git clone --single-branch -b release_5.0.0 https://github.com/ablab/quast.git
RUN cd /quast; ./setup.py install_full
ENV PATH $PATH:/quast

RUN conda install -c bioconda samtools
