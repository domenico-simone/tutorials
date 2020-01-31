FROM continuumio/miniconda3:4.7.10

RUN conda config --env --add channels conda-forge
RUN conda config --env --add channels bioconda
RUN conda create -n tutorials

# Activate environment
ENV PATH /opt/conda/envs/tutorials/bin:$PATH
ENV CONDA_DEFAULT_ENV tutorials
ENV CONDA_PREFIX /opt/conda/envs/tutorials

RUN echo "conda activate tutorials" >> ~/.bashrc
SHELL ["/bin/bash", "-c"]

# install tools
RUN conda install -c bioconda -c faircloth-lab fastqc sickle-trim scythe megahit prokka multiqc quast bowtie2 pilon samtools roary fasttree

# fancy CLI prompt :)
RUN echo "export PS1=\"\[\e[0m\e[47m\e[1;30m\] :: tutorials :: \[\e[0m\e[0m \[\e[1;34m\]\]\w\[\e[m\] \[\e[1;32m\]>>>\[\e[m\] \[\e[0m\]\"" >> /home/$CONTAINER_USER/.bash_profile

CMD /bin/bash -l