#!/bin/bash

tee /virus-discovery/config.json > /dev/null <<EOT
{
  "app": {
    "host": "0.0.0.0",
    "port": 5000,
    "output_path": "${OUTPUT_PATH}",
    "uploads_path": "${UPLOADS_PATH}",
    "uploads_size": ${UPLOADS_SIZE}
  },
  "db": {
    "host": "${DB_HOST}",
    "username": "${DB_USERNAME}",
    "password": "${DB_PASSWORD}",
    "schema": "${DB_NAME}"
  }
}
EOT

export PYTHONPATH=/usr/lib/python3/dist-packages/

python -u /virus-discovery/main.py