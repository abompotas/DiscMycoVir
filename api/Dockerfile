FROM python:3.9-bullseye

ENV OUTPUT_PATH=""
ENV UPLOADS_PATH=""
ENV UPLOADS_SIZE=8192
ENV DB_HOST=""
ENV DB_USERNAME=""
ENV DB_PASSWORD=""
ENV DB_NAME=""
ENV ADAPTER=""
ENV SLIDING_WINDOW=""
ENV MIN_LENGTH=""

WORKDIR /virus-discovery
COPY . /virus-discovery

RUN pip install -r requirements.txt && chmod 755 /virus-discovery/run.sh

CMD /virus-discovery/run.sh