# base_sequence_analysis

## Author
Byunghyun Ban

halfbottle@sangsang.farm

CTO / CBO of Imagination Garden Inc.


## Usage
Analysis on DNA base sequence to find transcriptable sites.

## Required Packages
None

## How to use

> import sequence_analizer as sa
> analizer = sa.Sequence_analizer("test.txt")

### find all possible genes and proteins
> analizer.find_all_possible_transcription()

### find all possible transcriptable sequences
> analizer.find_genes()

### find all possible proteins
> analizer.find_all_possible_proteins(

### Draw the gene and protein information
> analizer.protein
> analizer.transcriptable_sequences

### find all and save result
> analizer.save_result_as_FASTA_format()

Then results are save as separate files.

## Citation
Please Contact
>halfbottle@sangsang.farm