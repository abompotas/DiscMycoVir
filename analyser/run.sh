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
    "outfolder": "${OUTPUT_PATH}",
    "db2bcompared": "${PDB_PATH}",
    "kflag": ${K_FLAG},
    "dist": ${DIST},
    "max_mem": ${UPLOADS_SIZE},
    "seq_type": "${SEQ_TYPE}",
    "extensive": ${EXTENSIVE},
    "cpus": ${CPUS}
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
python -u /virus-discovery/process_queue.py