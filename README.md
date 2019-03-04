# dudeML
detection of duplications and deletions using machine learning. Note this tool works on most read data mapped to reference genome e.g. single-end short reads or MinION data, though all examples provided here use 100bp paired end reads.

# 1. Requirements

## Python-based:
* Python3
* pandas
* numpy
* scikit-learn
## External
* bedtools
* short-read simulator (e.g. wgsim)
* short-read aligner (e.g. BWA)
* short-read parser (e.g. SAMtools)
## Modifications
Within the script, the path for bedtools and wgsim need to be set. If both of these tools are within the path, they can be left as they are.

# 2. Functions
1. **winStat**
Finds the average coverage of a window in the chromosome, relative to the average chromosome coverage.
2. **fvecSample**
Reformats the bedfile with coverage information into sets of windows surrounding a focal window.
3. **fvecTrain**
Reformats the bedfile with coverage information into sets of windows surrounding a focal window. Also includes information on if a CNV is present in the window, and the estimated number of copies of that window per chromosome.
4. **simCNV**
Generate coordinates for random CNVs in the fasta file input, after accounting for repetitive content.
5. **recreateTotal**
If already known deletions and duplications are being used in the training data, this function skips the simulation of CNVs and instead generates a file with positions where CNVs should be, for simChrs.
6. **simChr**
Masks known repetitive content in the given chromosome and generates chromosomes with simulated CNVs.
7. **classify**
Given a training file (generated in formatTrain) and a sample file (generated in formatSample), will predicted windows with CNVs based on coverage and standard deviation of coverage.
8. **predict**
Given a training file (generated in formatTrain) and a sample file (generated in formatSample), will predicted windows with CNVs based on coverage and standard deviation of coverage.
9. **simReads**
Simulates read pairs of chosen length to a certain coverage of the chosen chromosome, requires WGsim.
10. **subTrain**
Downsample a training file by a certain percentage or to a certain number of each category. 
11. **summarize**
Combines called CNVs, and if known CNVs are provided, tells you if called CNVs are True-positives or otherwise.
12. **winStatExtra**
Creates summary windows based on known coverage estimates.
13. **covSummary**
Summarizes the coverage of each chromosome in the genomeCoverageBed file.

# 3. Input file formats
## Fasta
The reference sequences for mapping and for generating training files, in the following format:
>\>Chr1

>aagagcctatatca

>\>Chr2

>aagagcctatatca


## Coverage file
Generated from a BAM file of reads mapped to the reference fasta file, using the command:

	GenomeCoverageBed -d -ibam BAM > BED
Formatted as chromosome	position	coverage:
>Chr1	1	0

>Chr1	2	1

>Chr1	3	1

>Chr1	4	2


## Duplications
A bed file of known (or simulated) duplications, with the number of copies per chromosome and the frequency of the duplication in the sampled data (e.g. the number of chromosomes with this duplication/ the total number of chromosomes):
>Chr1	1000	1344	dup	3	1.0

>Chr1	2455	6700	dup	2	0.5

>Chr1	34501	36119	dup	2	1.0

>Chr1	45117	48932	dup	4	0.5

## Deletions
A bed file of known (or simulated) deletions, with the number of copies per chromosome and the frequency of the deletion in the sampled data (e.g. the number of chromosomes with this deletion/ the total number of chromosomes):
>Chr1	1000	1344	del	0	1.0

>Chr1	2455	6700	del	0	0.5

>Chr1	34501	36119	del	0	1.0

>Chr1	45117	48932	del	0	0.5

## Fvec files
A type of bed file, containing information on CNVs and copy number if a training file, and containing strain ID if a sample file. Reformats the bedfile from winStat to a feature vector, summarizing the windows around the focal window.
>2L	22240500	22240550	N	1.0	0.64	0.151	0.92	0.071	1.04	0.134	1.04	0.101	1.2	0.075	1.12	0.112	1.44	0.132	1.12	0.168	1.12	0.124	1.6	0.163	1.42	0.145

# 4. A simple walkthrough
## A. Simulate training data

We first downloaded the melanogaster reference genome and masked any repeats on the chromosome, before extracting chromosome 2L to use as a base for a training set. Following that, we simulated CNVs for a homozygous individual, requiring 1 set of chromosomes to be generated. Its important that the training data is as similar as possible to the sample being tested, so attempt to generate a file with similar coverage as your sample with a similar number of chromosomes (e.g. 2 for a heterozygote). The only files that need keeping after the chromosomes are simulated are the total files and del.1.bed and dup.1.bed, which contains the information about the simulated CNVs.

    repeatmasker -pa 4 -gff -gccalc -s -lib repbase.fa fasta/Dmel_iso1.fa
    bwa index fasta/Dmel_iso1.fa.masked
    
    for i in train test
    do
    mkdir ${i}_het
    python3 scripts/dudeML.py simCNV -fasta Dmel_iso1.fa.masked -cnvCount 10000 -d ${i}_sim -N 1
    python3 scripts/dudeML.py simChr -fasta Dmel_iso1.fa.masked -TE Dmel_iso1.gff -cnvBed ${i}_sim/total.1.bed -d ${i}_sim -id 1
    done


## B. Estimating coverage in training and test data
  We next simulated reads using WGSIM within dudeML and following mapping, used bedtools to calculate coverage per site, we then used a custom python script to find the mean coverage of each chromosome to find the relative coverages of each window. In this case a homozygote was simulated.
  
    for i in train test
    do
    python3 scripts/dudeML.py simReads -fasta Dmel_iso1.fa.masked -cov 20 -d ${i}_sim -RL 100 -id ${j}
    bwa mem -t 4 Dmel_iso1.fa.masked ${i}_sim/1_20_1.fq ${i}_sim/1_20_2.fq | samtools view -Shb - | samtools sort - > ${i}_sim/total.bam
    genomeCoverageBed -d -ibam ${i}_sim/total.bam > ${i}_sim/total.bed
    python scripts/dudeML.py winStat -i${i}_sim/total.bed -o ${i}_sim/total_50.bed -w 50 -s 50
    done

## C. Reformatting sample and training datasets.    
  We then removed repetitive regions, reformatted the data to show the relative coverage of the focal window and the 5 windows on each side. We also prepared the data to filter and extract the regions with known duplications or deletions in training the file. We also labelled CNVs in the test dataset for comparison later.

	python3 scripts/dudeML.py fvecTrain -i train_sim/total_50.bed -o train_sim/total_50train.bed -w 50 -TE Dmel_iso1.gff -dups train_sim/dup.1.bed -dels train_sim/del.1.bed  -windows 5
	python3 scripts/dudeML.py fvecSample -i test_sim/total_50.bed -w 50 -o test_sim/total_50sample.bed -id test_sim-TE Dmel_iso1.gff -windows 5

## D. Predicting CNVs using the generated files.  
Following this, you can create a classifier from one of the training features vector files generated and test out predictions of CNVs in the other file.

	python3 scripts/dudeML.py classify -i test_het/train_het/total_50train.bed -o train_het/total_50train.sav
	python3 scripts/dudeML.py predict -i test_het/total_50sample.bed -t train_het/total_50train.sav -o test_het/total_50pred.bed -windows 5

Alternatively, if multiple training files have been generated, these can be used to bootstrap the predicted CNVs, allowing you to take a consensus estimation of CNVs (much more conservative). In this case, the training file is set as a directory containing the training files.

    for i in {0..99}
    do
    python3 scripts/dudeML.py classify -i training/train_${i}.bed -o training/train_${i}.sav
    done
    
	python3 scripts/dudeML.py predict -i test_het/total_50sample.bed -t training/ -o test_het/total_50pred_bootstrap.bed -windows 5

## E. Using real data in dudeML
Real data with known structural variants can be used as a training set using the following pipeline.

	genomeCoverageBed -d -ibam knownCNV.bam > knownCNV.bed
	
	python dudeML.py winStat -i knownCNV.bed -o knownCNV_50.bed -w 50
	python3 dudeML.py fvecTrain -i knownCNV_100.bed -o knownCNV_50.bed -w 50 -TE repeats.gff -dups knownDUP.bed -dels knownDEL.bed

Real data can also be used as the test sample to identify unknown CNVs.

	genomeCoverageBed -d -ibam unknownCNV.bam > unknownCNV.bed
	
	python dudeML.py winStat -i unknownCNV.bed -o unknownCNV_50.bed -w 50
	python3 dudeML.py fvecSample -i unknownCNV_50.bed -o unknownCNV_50_sample.bed -w 50 -TE repeats.gff

The deletions and duplications bedfile should have 6 columns, showing the chromosome, start, end, type of CNV, frequency and number of copies of the CNV.
	
	2L  28779  30880  dup  0.5 2
	2L  41020  42111  dup  1.0 2
	2L  42277  42668  dup  0.5 2
	2L  55715  55717  dup  1.0 9
	2L  61325  61881  dup  1.0 2
	2L  69942  70335  dup  1.0 6
	2L  70571  72017  dup  1.0 5

Then the generated training and test sets can be used to find CNVs.
	
	python3 dudeML.py predict -i unknownCNV_50_sample.bed -t knownCNV_50_train.bed -o unknownCNV_50_pred.bed
