FROM ubuntu:22.04

ENV DISCVIR_WEBAPP="localhost:8100/#"
ENV DISCVIR_API="localhost"
ENV BATCH_SIZE=10
ENV SLEEP_SEC=20
ENV OUTPUT_PATH="/virus-discovery/output"
#----> PDB
ENV PDB_PATH=""

ENV K_FLAG=0
ENV DIST=10
ENV SASA_THRESHOLD=0.75
ENV DOCK_THRESHOLD=0.1
ENV EXTENSIVE=false
ENV CPUS=1
ENV DB_HOST="localhost"
ENV DB_USERNAME="root"
ENV DB_PASSWORD=""
ENV DB_NAME="virus_discovery_api"
ENV SMTP_SERVER="localhost"
ENV SMTP_PORT=465
ENV SMTP_ADDRESS=""
ENV SMTP_USERNAME=""
ENV SMTP_PASSWORD=""

WORKDIR /virus-discovery

ENV TZ="Europe/Athens"
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update \
&& apt install -y \
    software-properties-common \
    build-essential \
    dirmngr \
    curl \
    wget \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgfortran5 \
    mpich \
    gromacs \
    libfontconfig1-dev \
    python3-dev \
    python3-pip \
    libopenmm7.7 \
    libopenmm-dev \
    libopenmm-plugins \
    python3-simtk \
    openbabel \
    libopenbabel-dev \
    python3-openbabel \
    pymol \
    python3-pymol \
    python3-lxml \
    python3-distutils \
    pdb2pqr \
&& update-alternatives --install /usr/bin/python python /usr/bin/python3 1

#RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
#&& add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
#&& apt update \
#&& apt install -y r-base r-base-dev \
#&& R -e "install.packages('ggplot2')" \
#&& R -e "install.packages('gplots')" \
#&& R -e "install.packages('dplyr')" \
#&& R -e "install.packages('stringr')" \
#&& R -e "install.packages('igraph')" \
#&& R -e "install.packages('BiocManager'); BiocManager::install(version = '3.17'); BiocManager::install('annotate'); BiocManager::install('topGO'); BiocManager::install('biomaRt'); BiocManager::install('org.Hs.eg.db'); BiocManager::install('Rgraphviz')" \
#&& R -e "install.packages('GOxploreR')"


COPY . /virus-discovery

RUN python3 -m pip install lib/pdbfixer \
&& python3 -m pip install openbabel-wheel \
&& python3 -m pip install --no-deps -r requirements.txt

ENV ASSAM_V2=/virtuous-pocketome/lib/assam_linux_v2
RUN cp /virtuous-pocketome/lib/PLANTS/* /usr/local/bin/ \
&& chmod 755 /usr/local/bin/PLANTS1.2_64bit /usr/local/bin/SPORES_64bit \
    $ASSAM_V2/*.sh $ASSAM_V2/bin/aamc_pat $ASSAM_V2/bin/aamc_ssid $ASSAM_V2/bin/aspsettings $ASSAM_V2/bin/mastermod \
    $ASSAM_V2/bin/pdbaamc $ASSAM_V2/bin/aspreg2f $ASSAM_V2/bin/aspsup $ASSAM_V2/bin/matmult1r \
    /virtuous-pocketome/run.sh

CMD /virus-discovery/run.sh