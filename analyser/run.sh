#!/bin/bash

tee /virus-discovery/config.json > /dev/null <<EOT
{
  "webbapp": "${DISCVIR_WEBAPP}",
  "api": "${DISCVIR_API}",
  "queue": {
    "batch": ${BATCH_SIZE},
    "sleep": ${SLEEP_SEC}
  },
  "args": {
    "uploads": "${UPLOADS_PATH}",
    "output": "${OUTPUT_PATH}",
    "adapter": "${ADAPTER}",
    "sliding_window": "${SLIDING_WINDOW}",
    "min_length": ${MIN_LENGTH},
    "threads": ${THREADS},
    "max_memory": "${MAX_MEMORY}"
  },
  "db": {
    "host": "${DB_HOST}",
    "username": "${DB_USERNAME}",
    "password": "${DB_PASSWORD}",
    "schema": "${DB_NAME}"
  },
  "smtp": {
    "server": "${SMTP_SERVER}",
    "port": ${SMTP_PORT},
    "address": "${SMTP_ADDRESS}",
    "username": "${SMTP_USERNAME}",
    "password": "${SMTP_PASSWORD}"
  }
}
EOT
python -u /virus-discovery/process_queue_analysis.py &
python -u /virus-discovery/process_queue_trimming.py &
python -u /virus-discovery/process_queue_discovery.py