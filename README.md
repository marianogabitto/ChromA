# ChromA: Chromatin Landscape Annotation Tool

ChromA is a probabilistic model to annotate chromatin regions into accessible or inaccessible, open or closed, based on their ATACseq profile. ChromA can process bulk datasets, single-cell or integrate information from a combination of both. Even more, ChromA can integrate information from different replicates or different cellular populations to create a consensus representation of chromatin accessibility.

At the moment, we are building a webpage with extensive documentation. We give a brief introduction here on how to install and run ChromA but contact the authors or open an issue in github if you are encountering errors.

Publication: https://www.nature.com/articles/s41467-020-14497-5

# Installation
The current version of ChromA is implemented in Python 3 and can be installed by running:

    pip install git+https://github.com/marianogabitto/ChromA@v2.1.2

ChromA requires the following packages. In principle, the pip installer should manage their installation but just in case, here is a command that will install all the dependencies.

    pip install matplotlib==3.0.2 nose==1.3.7 psutil pysam==0.15.2 ray==0.6.3 scipy==1.2.0 seaborn==0.9.0 setproctitle==1.1.10 numpy==1.15.4

You can verify the correct installation of ChromA by running:

    ChromA

The output should look like this:
    ChromA
>usage: ChromA [-h] [-v] [-ve VERBOSE_LEVEL] [-sp SAVE_POST]
              {dnase,atac,Cutrun,count,filter} ...

>ChromA

>positional arguments:
  {dnase,atac,Cutrun,count,filter}
                        help for sub-command
    dnase               Run ChromA for DNASE experiments.
    atac                Run ChromA for ATAC experiments.
    Cutrun              Run ChromA for CutRun experiments.
    count               Count Fragments in peaks.
    filter              Filter Fragments given barcodes.

>optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -ve VERBOSE_LEVEL, --verbose VERBOSE_LEVEL
                        verbose level. 1 display messages, 0 omit them
                        (default is 0).
  -sp SAVE_POST, --posterior SAVE_POST
                        Save Posterior State. This is memory intensive.
                        (default=False, Not Required).

>ChromA: error: the following arguments are required: -i/--input, -sb/--saveBedFile

Omit the error line as it is telling you that no input file was detected. However, a potential error message can be displayed "Failed to Load Cpp Core" . If this is the case, please refer to troubleshooting ChromA Installation.


# Get Started:
ChromA is an easy to use software tool to create annotations of chromatin landscape from ATAC-seq information. In what follows, we detail the necessary steps in order to run ChromA on your datasets.

CHECK CHROMA INSTALLATION:  Let's first verify ChromA's correct installation. In the next lines, we describe commands to be run with >, and desired output with italics.

    ChromA -v
> ChromA 2.1.0

RUNNING CHROMA ON A RESTRICTED CHROMOSOMAL SUBSET:  Next, verify that your file meets the FILE REQUIREMENTS. Once you have proper files, you can run ChromA using:

    ChromA atac  -spec mm10 -i "my.sorted.index.bam" -sb output -reg True

RUNNING CHROMA ON THE ENTIRE GENOME: To run ChromA on the entire genome please run:

    ChromA atac  -spec mm10 -i "my.sorted.index.bam" -sb output

We support the following genomes: mm10 / hg38 / hg19/ dm6. In case you want a new genome being supported, please contact the authors.

# File Requirements
ChromA supports different file formats:
* BAM files. BAM files should be indexed and sorted before running ChromA. The corresponding index file (.bai) should be located on the same directory where the bam file is present.
* TSV files. ChromA supports tab separated files. The files have no header and the four columns are: chrs, start, stop, #duplicates.
* TSV.GZ files. ChromA supports tab index files (tabix) such as the output of cell ranger. This file tipycally contains the out put of a single-cell experiment. As such, it is a compressed file resulting from a 4-column tab separated no-header file. The four columns are chr, start, stop, #duplicates.  

# Notes
If running ChromA returns an error in which "libfwdbwdcpp.so" is mentioned, this is due to troubles finding the C++ library that performs calculations. Please refer to our troubleshoot section to overcome this issue.
If running ChromA returns an error in which "worker.py" is mentioned, this could be caused by our parallelization support. Please refer to our troubleshoot section to overcome this issue.

RUNNING CONSENSUS CHROMA: To run Consensus ChromA integrating information from different replicates, integrate information from bulk and single cell or to create a common representation between different populations, please list all the files after the -i command. Be aware that this command will take longer to run.

    ChromA atac -i my.sorted1.bam my.sorted2.bam my.tn5binding.tsv.gz -reg True  -sb output.bed -spec mm10

RUNNING CHROMA ON A CLUSTER USING SLURM:To run ChromA on a computer cluster using the slurm scheduler, we create a bash script named run.sh. run.sh has the following structure:

    cat run.sh
    
> #!/bin/bash

The following lines load cuda support and select the python 3 environment in which ChromA is installed

> module load gcc/7.3.0

> module load cuda/8.0.61

> module load cudnn/v6.0-cuda-8.0

> source /mnt/home/mgabitto/python_3/bin/activate

The following lines initialize ray without gpus and then run ChromA

>CUDA_VISIBLE_DEVICES=''

> ChromA atac  -spec mm10 -i "my.sorted.index.bam" -sb output

To run our bash script on the slurm scheduler please run:

    sbatch -N 1 --exclusive -t 0-10:00:00 run.sh

# CHROMA OUTPUTS:
ChromA generates 4 files as its output. In general, ChromA minimizes output to the console, logging all the information to different files. After a successful run, you should find the following 4 files:

    ChromA atac  -spec mm10 -i "my.sorted.index.bam" -sb output

Within the directory where my.sorted.bam resides:
* my.sorted_insert_size.pdf

Within the directory where the ChromA command is invoked:
* output_allpeaks.bed
* output.log
* output_metrics.log

my.sorted_insert_size.pdf creates an insert size distribution graphs for dataset quality control.
output_allpeaks.bed is a bed file that contains accessible regions.
output.log logs all the information of the current run.
output_metrics.log generates dataset metrics (such as FRIP, SN, #Reads) for quality control.
