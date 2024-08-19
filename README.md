# Virus Discovery


## Requirements

- git
- Docker


## Quick Run

### Installation

If this is the first run, follow the steps below before starting the app.

1. Download the code:
    ```shell
    git clone https://gitlab.com/ceid-makris-team/bioinformatics/virus-discovery
    ```

2. Edit the file ```example.env``` that contains the configuration options and rename it to ```.env```:
    ```shell
    nano example.env
    mv example.env .env
    ```

3. Create the persistent storage volumes:
    ```shell
    docker volume create virus-discovery-db
    docker volume create virus-discovery-outputs
    docker volume create virus-discovery-uploads
    ```
   (Optional) For performance reasons it is advised to run the app in machine with SSD storage. 
   However, these volumes can grow in size pretty quickly and are not critical for the app's performance.
   If you wish to create them in a separate storage (e.g. HDD or NAS) you can use the 
   __[--opt](https://docs.docker.com/reference/cli/docker/volume/create/#opt)__ 
   parameter of the ```docker volume create``` command:
    ```shell
    docker volume create --opt type=none --opt o=bind --opt device=/hdd/path/to/folder virus-discovery-db
    docker volume create --opt type=none --opt o=bind --opt device=/hdd/path/to/folder virus-discovery-outputs
    docker volume create --opt type=none --opt o=bind --opt device=/hdd/path/to/folder virus-discovery-uploads
    ```


### Run

Run the application
```shell
docker compose up -d
```
Subsequently, if you haven't altered the URL settings in the ```.env``` you can access the app through 
[http://localhost:8000](http://localhost:8000) (or whatever you have specified)

### Stop

Stop the application
```shell
docker compose down
```