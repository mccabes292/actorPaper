# Code for obtaining the real dataset

* `getSamp1.sh` - Code for pulling the FASTQ files from the list of IDs.
* `salmon1.sh` - Run `salmon` for transcript quantification. Needs to be run per ID.
* `quantToCounts.R` - Combines output from `salmon` into one matrix containing all samples.
