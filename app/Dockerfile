FROM nginx

ENV DISCVIR_API="http://localhost:5000"

WORKDIR /virus-discovery
COPY ./www /virus-discovery

COPY ./run.sh /run.sh
RUN chmod +x /run.sh

CMD ["/run.sh"]
