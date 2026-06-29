## fq2vcf

This is a rewrite of the mapping/calling pipeline, done originally on CSC Mahti. With minor modifications it should run on other environments.

### Running the script 

An example of running the script is here: 
```
#!/bin/bash
#SBATCH --job-name=run_nf
#SBATCH --account=project_200xxxx
#SBATCH --time=72:00:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=34
#SBATCH --output=slurm/run_nf-%j.out

set -euo

export PATH=/projappl/project_200xxxx/bin:$PATH
module load nextflow

export NF_AWS_ACCESS_KEY="$(grep access_key ~/.s3cfg_200xxxx | awk '{print $3}' | tr -d ' ')"
export NF_AWS_SECRET_KEY="$(grep secret_key ~/.s3cfg_200xxxx | awk '{print $3}' | tr -d ' ')"
export TOWER_ACCESS_TOKEN="$(cat ~/.token)"
export NXF_VER=25.10.2

nextflow run map_and_call.nf --samples three_samples.yaml -with-tower
```

The script doesn't load or provide paths to programs and they need to be made available by the user (see `PATH` line above).

The script is independent of `s3cmd` and reads/writes Allas natively. It gets the S3 parameters from a config file (above, `~/.s3cfg_200xxxx`). The S3 parameters should be for the Allas project that hosts the data or has permissions to the input data and the output data bucket. The project doesn't need to be the active project. If using `-with-tower` (see below), the token needs to available in a file or integrated in the script.

### Defining samples with FASTQ data 

The sample data is given in a YAML file. This file defines two samples with pair-end short-read data (`type: SR`):
```
samples:
  - sample_id: PV-K40
    fastq:
      - s3://fastq_data/harborseal/PV-K40/PV-K40_1.fq.gz
      - s3://fastq_data/harborseal/PV-K40/PV-K40_2.fq.gz
    bam:  s3://processed_ribbon/harborseal/PV-K40/PV-K40_markdup.cram
    gvcf: s3://processed_ribbon/harborseal/PV-K40/PV-K40.gvcf.gz
    type: SR
  - sample_id: PV-K48
    fastq:
      - s3://fastq_data/harborseal/PV-K48/V350160768_L01_507_1.fq.gz
      - s3://fastq_data/harborseal/PV-K48/V350160768_L01_507_2.fq.gz
      - s3://fastq_data/harborseal/PV-K48/V350160768_L02_507_1.fq.gz
      - s3://fastq_data/harborseal/PV-K48/V350160768_L02_507_2.fq.gz
      - s3://fastq_data/harborseal/PV-K48/V350160768_L03_507_1.fq.gz
      - s3://fastq_data/harborseal/PV-K48/V350160768_L03_507_2.fq.gz
      - s3://fastq_data/harborseal/PV-K48/V350160768_L04_507_1.fq.gz
      - s3://fastq_data/harborseal/PV-K48/V350160768_L04_507_2.fq.gz
    bam:  s3://processed_ribbon/harborseal/PV-K48/PV-K48_markdup.cram
    gvcf: s3://processed_ribbon/harborseal/PV-K48/PV-K48.gvcf.gz
    type: SR
```
The accepted suffixes for pair matching are listed in `preprocess.nf/ALIGN_BWA`.

This defines one HiFi-sample (`type: HIFI`):
```
samples:
  - sample_id: Tuukka
    fastq:
      - s3://fastq_DNA/pacbio/Saimaa/Phs499_Tuukka/HiFi_bam/m64145_231127_092254.bcAd1043T--bcAd1043T.bam
      - s3://fastq_DNA/pacbio/Saimaa/Phs499_Tuukka/HiFi_bam/m64145_231127_092254.bcAd1053T--bcAd1053T.bam
      - s3://fastq_DNA/pacbio/Saimaa/Phs499_Tuukka/HiFi_bam/m64145_231224_080427.bcAd1053T--bcAd1053T.bam
      - s3://fastq_DNA/pacbio/Saimaa/Phs499_Tuukka/HiFi_data/m84212_240209_092418_s1.hifi_reads.bcAd1053T.bam
    bam:  s3://processed_ribbon/Saimaa_ringedseal/Tuukka/Tuukka_hifi.cram 
    gvcf: s3://processed_ribbon/Saimaa_ringedseal/Tuukka/Tuukka_dv.gvcf.gz 
    type: HIFI
```
HiFi data can be either as FASTQ(.gz) or BAM, but the two formats cannot be mixed in one sample.

### Defining samples with CRAM data 

If the mappping has worked but the calling not, the script should be able to continue from a CRAM file in Allas. The file is defined in the field `bam:`. The field `fastq:` has to be defined but no files listed below it.

 This file defines two samples with CRAM files consisting of pair-end short-read data (`type: SR`):
```
samples:
  - sample_id: PV-K40
    fastq:
    bam:  s3://processed_ribbon/harborseal/PV-K40/PV-K40_markdup.cram
    gvcf: s3://processed_ribbon/harborseal/PV-K40/PV-K40.gvcf.gz
    type: SR
  - sample_id: PV-K48
    fastq:
    bam:  s3://processed_ribbon/harborseal/PV-K48/PV-K48_markdup.cram
    gvcf: s3://processed_ribbon/harborseal/PV-K48/PV-K48.gvcf.gz
    type: SR
```

The functionality is not defined for HiFi data but it would be fairly simple to do.


### Analysis programs included in the script 

The programs run can be seen in the code and can be guessed from the process names. The SR-mode runs processes:
`BWA_INDEX; ALIGN_BWA; REALIGN_INDELS; MERGE_BAMS; RUN_HAPLOTYPECALLER; MERGE_GVCFS; CONVERT_CRAM_AND_PUBLISH; PUBLISH_GVCF`
and the HiFI-mode runs processes:
`ALIGN_MINIMAP2; RUN_DEEPVARIANT; CONVERT_CRAM_AND_PUBLISH; PUBLISH_GVCF_VCF`

The BWA index is published and stored in `reference` folder. It needs to be done only once.

### Settings for the reference genome

Settings specific to the reference genome are in the beginning of the `main.nf` file.

```
params.reference       = 's3://processed_ribbon/reference/Ribbon_PR_18_6_2024.fa.gz'
params.refsubsets      = 's3://processed_ribbon/reference/ref_subsets.tgz'
params.cramreference   = '/scratch/project_200xxxx/references/Ribbon_PR_18_6_2024.fa.gz'
```

The genome from `params.reference` is copied to `params.cramreference` and hardcoded in the CRAM files. (This is useful but not too serious. If the hardcoded reference is not there, CRAM files can be opened by providing the reference in the `samtools` command.)

The `params.refsubsets` contains the subsets of the genome used for the parallelisation of otherwise unparallelisable processes. The subset files should be named as `S0?.list`-`S??.list` so that they get naturally sorted. The file above contains subsets S01-S10, but a smaller number of subsets should work. Ideally, the subsets are of similar size and take roughly the same time to run.

The script will not run without the subsets. If subsetting is not wanted, one can create a fake subset listing all chroms/contigs: `cut -f1 reference.fa.fai > S01.list && tar czf fake_subset.tgz S01.list`.   


### CPU and RAM allocation

The CPU and RAM requirements of the processes are listed in `nextflow.config`. On Mahti, RAM cannot be allocated and each core gives 1.875GB of RAM. Above, 34 cores is enough to run BWA for 3 parallel samples, each requiring 10 cores and 20GB of RAM. (Technically, 32 cores should be enough: `32*1.875GB=60GB`, but the script needs some memory for other tasks as well.) However, it cannot run 3x10 parallel jobs of `REALIGN_INDELS` or `RUN_HAPLOTYPECALLER` as those require 2 cores and 4GB each. (That would require `3x10x4GB=120GB` while we only have `34*1.875GB=63.75GB`.)   

Of the processes, `REALIGN_INDELS` and `RUN_HAPLOTYPECALLER` are parallelised with the genome subsets, the number of subsets defined by the user. All other processes are run for the full data per sample. Look at the CPU and RAM usage of each process in `nextflow.config` and calculate your allocation accordingly. Note that, on Mahti, the RAM limit is often more important than the core limit! 

### Nextflow Tower

To monitor the jobs through the Tower interface, one needs to create a user account at https://cloud.tower.nf/ and then a token. The token can be created through the account icon in the top-right corner.


### Creation of sample files

A simple script is provided for the creation of sample YAML files. It may not work for all cases but gives an idea how the file is structured.

```
cat > three_samples.txt << EOF
PV-K40
PV-K48
PV-K49
EOF
```

```
bash create_sample_file.sh three_samples.txt harborseal > three_samples.yaml
```

Another script is available for cases where the file paths are given in a local file.

## Jointcalling

An attempt has been started for the join-call of GVCF data. 
