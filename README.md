# NGS Toolkit

This repository is developed for analyzing NGS datasets, specifically the ones related with genome-wide DNA Damage and Repair.

## Setup
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
docker-compose run main bash
```

`docker-compose` maps the data directory to `/data` in the container. All the programs should be installed for the docker container.

Now it is time to build a reference genome index for the alignment program. Here we use `Bowtie2`.

```
/NGStoolkit/Docker/prepareReferenceGenome.sh
```
This will take some time, feel free to have a cup of coffee.

After the reference genome index is built successfully, we can run our pipeline.

```
cd data
/NGStoolkit/stable/
/NGStoolkit/stable/XR-seq-basics.sh
```

Here we go!

## Authors
  * Ogun Adebali

## Licence
  This project is licensed under the MIT License
