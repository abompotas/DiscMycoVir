# Virus Discovery

## Requirements

- git
- Docker


## Quick Run

### Installation

If this is the first run, follow the steps below before starting the app.

Download the code:
```shell
git clone https://gitlab.com/ceid-makris-team/bioinformatics/virus-discovery
```

Edit the file ```example.env``` that contains the configuration values and rename it to ```.env```:
```shell
nano example.env
mv example.env .env
```

Create the persistent storage volumes:
```shell
docker volume create virus-discovery-db
docker volume create virus-discovery-outputs
docker volume create virus-discovery-uploads
```

### Run

Run the application
```shell
docker compose up -d
```

### Stop

Stop the application
```shell
docker compose down
```