FROM python:3.9-bullseye

ENV UPLOADS_PATH="/virtuous-pocketome/uploads"
ENV UPLOADS_SIZE=64
ENV DB_PASSWORD=""

WORKDIR /virtuous-pocketome
COPY . /virtuous-pocketome

RUN pip install -r requirements.txt && chmod 755 /virtuous-pocketome/run.sh

CMD /virtuous-pocketome/run.sh