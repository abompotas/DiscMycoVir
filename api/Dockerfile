FROM python:3.9-bullseye

ENV UPLOADS_PATH=""
ENV UPLOADS_SIZE=8192
ENV DB_PASSWORD=""

WORKDIR /virus-discovery
COPY . /virus-discovery

RUN pip install -r requirements.txt && chmod 755 /virus-discovery/run.sh

CMD /virus-discovery/run.sh