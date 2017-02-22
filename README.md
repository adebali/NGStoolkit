# Sancar Lab Scripts

This repository is developed for project-specific NGS analysis.

## Prerequisites
  * bowtie
  * ...

## Setup
  ```bash
  cd ~
  git clone https://github.com/adebali/sancarLabUtils.git
  ./setup.sh
  ```
## Examples

## Running the tests

## Built with
  * python
  * LSF

## Authors
  * Ogun Adebali

## Licence
  This project is licensed under the MIT License

## Acknowledgments
  * ...


## Examples
### DamageSeq MiSeq
```
cd /proj/sancarlb/users/ogun/scripts/projects/DamageSeq/miseq
loopThroughFiles.py -code "bsub python pipeline.py #IN" -files "dataDir/160831_UNC23_0048_000000000-ARR90/*.fastq"

```

