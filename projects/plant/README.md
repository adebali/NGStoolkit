# Genome-wide Excision Repair Maps of _Arabidopsis thaliana_
Q
We generated XR-seq data sets collected at different circadian time points. Here is the pipeline package we generated for the analysis.

Reference genome: TAIR10

## External dependencies of the pipeline
* cutadapt
* bowtie
* bedtools
* pybedtools (python module)


## Usage

### Running determined pipeline parts

First, make sure the runFlags in the branches and each execution command are set based on your design. Then apply the following:

```python pipeline.py run -n 1```

1 is an example sample id in the ```samples.csv```

If you are on a cluster having ```slurm``` as a job scheduler you can run the array job submission:

```python pipeline_array -g 0```

-g option represents the group(s) to run (see ```samples.csv```)



### Running everything from scratch (rarely used)

```python pipeline.py run -n 1 --fromScratch```

