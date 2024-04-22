#!/bin/bash

tee /virtuous-pocketome/config.json > /dev/null <<EOT
{
  "webbapp": "${POCKETOME_WEBAPP}",
  "api": "${POCKETOME_API}",
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
    "sasaThreshold": ${SASA_THRESHOLD},
    "dockThreshold": ${DOCK_THRESHOLD},
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

python -u /virtuous-pocketome/process_queue.py