# NGS Toolkit

This repository is developed for analyzing NGS datasets, specifically the ones related with genome-wide DNA Damage and Repair.

## Setup

Custom setup assumes that you have the necessary programs installed and they are exectuble in your `$PATH`. 

The list of necessary programs are:
* Bowtie2
* cutadapt
* sra-tools
* samtools
* bedtools
* bedGraphToBigWig (uscstools)

To place the source code of this repository in your path, please follow these commands:
```
cd ~
git clone git@bitbucket.org:adebali/NGStoolkit.git
cd NGStoolkit
bash setup.sh
```

## Docker setup
Make sure that you have `Docker` and `docker-compose` are installed.

```
docker-compose -p ngs up --build -d
```

### Work in the container
`docker-compose` maps the data directory to `/data` in the container. All the programs should be installed for the docker container.

Go into the container with `docker-compose run main bash`

Now it is time to build a reference genome index for the alignment program. Here we use `Bowtie2`. Download and prepare the reference genome with `/NGStoolkit/Docker/prepareReferenceGenome.sh`. This will take some time, feel free to have a cup of coffee.

After the reference genome index is built successfully, we can run our pipeline `cd data && /NGStoolkit/stable/XR-seq-basics.sh`. Here we go!

### Work outside of the container

* Download and prepare the reference genome with `docker-compose run main /NGStoolkit/Docker/prepareReferenceGenome.sh`
* Run the pipeline with `docker-compose run main /NGStoolkit/stable/XR-seq-basics.sh`

## Work with your own data

* Move your `.fastq` file into the `data` directory.
* Edit the `XR-seq-basic.sh` and repalce the `SAMPLE` variable with the base sample name in your file. For example if you file is named as `myFile.fastq` the base name will be `myFile`. 
* If you want to retrieve the existing data set from SRA please see the `fastq-dump` command and replace the SRA acccession number with the one of interest. If you use youw own file please comment out that two lines in `XR-seq-basics.sh`.

## Authors
  * Ogun Adebali

## Licence
  This project is licensed under the MIT License
