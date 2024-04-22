# VirtuousPocketome

This tool was developed within the Virtuous Project (https://virtuoush2020.com/)

[![Virtuous button][Virtuous_image]][Virtuous link]

[Virtuous_image]: https://virtuoush2020.com/wp-content/uploads/2021/02/V_logo_h.png
[Virtuous link]: https://virtuoush2020.com/


### Repo Structure
The repository is organized in the following folders:

- code/
>Collecting the VirtuousPocketome core codes and programs

- data/
>Collecting example configuration files to run the VirtuousPocketome tool

- notebook/
>Collecting jupyter notebooks used to prepare the analysis and the figure added in the relative publication


### Authors
1. [Lorenzo Pallante](https://github.com/lorenzopallante)
2. [Marco Cannariato](https://github.com/marcocannariato)
3. [Lampros Androutsos](https://github.com/lamprosandroutsos)
4. [Eric Zizzi](https://github.com/ericzizzi) 


----------------
## Prerequisites
----------------

1. Create conda environment:

        conda create -n cbp python=3
        conda activate cbp

2. Install required packages:

        conda config --add channels conda-forge
        conda config --set channel_priority strict
        conda install matplotlib mdanalysis plip pdbfixer requests gromacs
        pip install pdb-tools art suds
    
    Note: in case you want to install jupyter, it is suggested to do it before the pip installation lines
    
3. Clone the `CBP` repository from GitHub

        git clone https://github.com/lorenzopallante/CBP

4. Install required additional softwares: 
   
   - GROMACS
    -> easy way: `sudo apt install gromacs` (Recommended way for the present application)
    -> hard way: compile it from source: https://manual.gromacs.org/documentation/current/install-guide/index.html
   
   - PLANTS and SPORES
    -> executable present in this repo 
    -> if you need them, download from http://www.tcd.uni-konstanz.de/

    - R (you should encounter some permission errors.. add permission to the R installation folders)
    `sudo apt-get update`
    `sudo apt-get install r-base`
    install additional packages from the R console (type 'R' in the terminal)
    `install.packages("ggplot2")` 
    `install.packages("gplots")` 
    `install.packages("dplyr")` 
    `install.packages("stringr")`
    `install.packages("GOxploreR")`
    `install.packages("igraph")`

    install topGO and relative additional packages using the following coding line: 
    `if (!require("BiocManager", quietly = TRUE))`
    `   install.packages("BiocManager")`
    `BiocManager::install(version = "3.17")`
    `BiocManager::install("topGO")`  
    `BiocManager::install("biomatRt")`
    `BiocManager::install("org.Hs.eg.db")`
    `BiocManager::install("Rgraphviz")`


    You may need to install these libs on your Lunux machine to be able to install BiocManager and related packages: 
    curl: `sudo apt-get install curl`
    libssl-dev: `sudo apt-get install libssl-dev`
    libcurl: `sudo apt-get install libcurl4-openssl-dev`
    xml2: `sudo apt-get install libxml2-dev`

5. Download the database folder collecting the PDB to be screened (approx. 80 GB)


Enjoy!        

-------------------------------
## How to use VirtuousPocketome
-------------------------------

The main code to run VirtuousPocketome is `master.py` within the code folder. 

To learn how to run, just type: 

    python master.py --help

The best way to run the code is to use a configuration file (example available in the data folder)

    python master.py --config config_template.txt

You need to provide at least the path to the protein-ligand complex pdb file and the chain labels identifying the two elements in the system. 

Briefly, the master.py code will call in order: 

1- `generate_assam_motifs.py` to generate the ASSAM motifs from the input file(s)
    
2 - `run_assam.py` to perform the Similarity Search using ASSAM

3 - `analyse_assam.py` to perform the multi-step filtering process 

The function library for the codes is `motif_search.py`.
    
------------------
## Acknowledgement
------------------

The present work has been developed as part of the VIRTUOUS project, funded by the European Union’s Horizon 2020 research and innovation program under the Marie Sklodowska-Curie-RISE Grant Agreement No 872181 (https://www.virtuoush2020.com/).
