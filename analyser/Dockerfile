FROM ubuntu:22.04

ENV DISCVIR_WEBAPP="localhost:8100/#"
ENV DISCVIR_API="localhost"
ENV UPLOADS_PATH="/virus-discovery/uploads"
ENV OUTPUT_PATH="/virus-discovery/output"
ENV SEQ_TYPE="fa"
ENV UPLOADS_SIZE=8192

ENV BATCH_SIZE=10
ENV SLEEP_SEC=20
ENV ADAPTER=""
ENV SLIDING_WINDOW=""
ENV MIN_LENGTH=""
ENV THREADS=1
ENV MAX_MEMORY=4096

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

RUN chmod 755 /virus-discovery/lib/virus_discovery_pipeline.sh /virus-discovery/run.sh

CMD /virus-discovery/run.sh