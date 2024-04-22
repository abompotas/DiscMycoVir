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
    curl \
    wget \
    libcurl4-openssl-dev \
    libssl-dev \
    python3-dev \
    python3-pip \
    fastqc \
    trimmomatic \
    trinityrnaseq \
    bwa \
    ncbi-blast+ \
&& update-alternatives --install /usr/bin/python python /usr/bin/python3 1

COPY . /virus-discovery

RUN python3 -m pip install -r requirements.txt

RUN chmod 755 /virtuous-pocketome/lib/virus_discovery_pipeline.sh /virtuous-pocketome/run.sh

CMD /virus-discovery/run.sh