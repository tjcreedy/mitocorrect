# Mitocorrect
Mitocorrect is a python program for correcting mitochondrial genome annotations.

## Installation
Download and run the scripts from this repository.

### Dependencies
You must have the biopython module installed, you can install this through pip, conda or your 
package manager.

You must also have [mafft](https://mafft.cbrc.jp/alignment/software/) installed and available on 
the PATH.

### System requirements
These scripts were developed and tested on Ubuntu Linux. 

## Preparing to run mitocorrect
### Mitogenome sequences
Mitocorrect works on complete or incomplete mitogenomes, with all or some genes, and can deal with
genes that are truncated by the end of a mitochondrial contig or are have unknown sections due to
poor sequencing coverage. Mitocorrect recognises truncations and does not apply length-based 
scoring to these annotations; however where internal sections of the gene are unknown, it is 
important to ensure that a roughly accurate number of Ns are used to ensure length-based scoring
is accurately applied.

### Annotating mitogenomes
Mitocorrect works on mitochondrial genomes that are already annotated. If you have unannotated 
mitogenomes, you can annotate them by using one of the many annotation tools out there. For small
quantities of mitogenomes, the [MITOS WebServer](http://mitos2.bioinf.uni-leipzig.de/index.py) is
pretty good. For larger quantities, you could use [MitoZ](https://github.com/linzhi2013/MitoZ/).

The annotated mitogenomes must be in 
[genbank format](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html). Mitocorrect can take as
input any number of genbank files, so your sequences can be in one or many files. Note that the 
genbank files are read using biopython so the format must conform to the biopython expectations. 
In some cases, it seems MitoZ omits the Division code and so does not form the LOCUS line 
correctly, which causes reading the genbank file to fail. In this case you will need to correct 
the genbank file before running mitocorrect.

### Preparing profiles
One way mitocorrect assesses potential annotations is through comparison with a set of reference
gene alignments, which we refer to as "profiles". These need to be constructed separately prior
to running mitocorrect, which expects a tab-delimited table giving the gene name and path for each
alignments. These alignments must comprise complete gene sequences from start to stop codons 
(including stops, these may be partial), and the species should be representative of the diversity
of input species. We recommend at least 30 sequences per gene, and we recommend performing amino
acid alignments as these are going to be more appropriate for identifying the positions of start
and stop codons.

### Specifications
Mitocorrect is designed to apply expert knowledge about annotation placement on a mitogenome to 
large quantities of mitogenome sequences. This expert knowledge is communicated to mitocorrect
in a specifications table, which describes the expectations mitocorrect should have about the 
position of a gene. Each gene should have two lines: one each for the start and stop. The table
should be a tab-delimited file with the following columns:

| Column             | Required | Default value | Notes                                                                                                                                                                                                                                            |
|--------------------|----------|---------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| gene               | yes      |               | The short name of the gene, should be repeated on every line                                                                                                                                                                                     |
| length             | yes      |               | The expected length of the gene, in nucleotide base pairs. Can be left blank if already specified for this gene on a previous line                                                                                                               |
| lengthvariation    | yes      |               | The expected percent variation of the gene from the supplied length, outside of which potential annotations will be scored down                                                                                                                  |
| end                | yes      |               | The end of the gene this row refers to, i.e. `start` or `stop`                                                                                                                                                                                   |
| overlap            | see note | *None*        | A `/`-separated list of 'context' genes (e.g. `NAD6,1/TRNP,-3) to use for assessing possible start/stop codons, and the number of bases of overlap expected from that gene. Positive values denote overlap, negative values denote gap.          |
| overlapmaxdistance | no       | `50`          | The maximum number of bases around the specified overlap position to search for possible start/stop codons                                                                                                                                       |
| searchcode         | no       | `N`           | Whether the list of start/stop codons to search for are `A`mino acids or `N`ucleotide sequences                                                                                                                                                  |
| searchsequence     | see note | *None*        | A `/`-separated list of sequences to search for (e.g. `ATG/ATT/ATA`) to find start/stop codons. These can be one or more bases, covering multiple codons: potential annotations will always cover the entire searchsequences. Can be amino acids |
| searchdistance     | no       | `30`          | The maximum number of bases around the current start/stop position of the annotation to search for possible start/stop codons                                                                                                                    |
| alignbody          | no       | `0`           | On the profile alignment, the number of bases from the modal start/stop to the start of the main body of the alignment. Positive values are in the 3' direction, negative values in the 5' direction. See explanation below                      |
| alignweight        | no       | `1.01`        | Value by which to weight the alignment score relative to other scores for the final scoring of potential annotations                                                                                                                             |
| positionweight     | no       | `1`           | Value by which to weight the position score relative to other scores for the final scoring of potential annotations                                                                                                                              |
| indelweight        | no       | `0.5`         | Value by which to weight the indel score relative to other scores for the final scoring of potential annotations                                                                                                                                 |

If a gene end does not have values for overlap and searchsequence, this will be skipped and the
existing start or stop codon will be retained.

We supply a minimal default version of this table in this repository. We expect that this will
not work optimally for most cases, but that it should be modified based on prior knowledge about
the taxon, and likely through iteratively running mitocorrect and inspecting the outputs. We 
mostly use mitocorrect on the beetles (Coleoptera), and provide our table for the beetles too. 
We would welcome receiving specification table from other taxa.

### TODO: write something about alignbody

## Running mitocorrect
### Taxonomic diversity
If you have highly diverse sequences (say, across multiple orders or classes), we recommend running
taxonomic groups (order, family) separately so as to ensure that the profiles used are 
representative.

### Multithreading
Mitocorrect can take advantage of multiple CPU cores/threads to correct mitogenomes in parallel.
This offers substantial performance benefits and is strongly recommended. Note that mitocorrect 
runs four output helper threads alongside one or more correction threads, so with `-t 1` you
should expect to see 5 mitocorrect processes - this is normal.

### Command reference
Mitocorrect is a simple python script, that can be run as follows to get usage information:

`python3 mitocorrect.py --help`

#### Required arguments

| Argument                 | Description                                                                                                                                                                                                                                               |
|--------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-s`/`--specifications`  | The path to a specifications file (see above) in tab-separated variables format                                                                                                                                                                           |
| `-g`/`--genbank`         | The path to one or more genbank-format files containing annotated mitogenomes to correct. Unless the argument `-1`/`--onefile` is used (see below), corrected sequences will be written to the output folder in the same set of files with the same names |
| `-a`/`--alignmentpaths`  | The path to a two column tab-separated variables file giving the name and path to a profile alignment (see above) for each gene                                                                                                                           |
| `-b`/`--translationtable`| The appropriate numeric genetic code translation table (see [here](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi))                                                                                                                             |

#### Common optional arguments

| Argument                 | Description                                                                                                                                                                                                               | Default                                             |
|--------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------|
| `-t`/`--threads`         | The number of simultaneous threads to use for correcting mitogenomes                                                                                                                                                      | `1`                                                 |
| `-c`/`--alignmenttype`   | The type of sequence data used in the supplied profile alignments (`aa` or `nt`)                                                                                                                                          | `nt`                                                |
| `-o`/`--outputdirectory` | The path to a directory (created if needed) to which intermediate and output files will be written. Note it is not recommended that this be set to the same directory holding the input files, as they may be overwritten | `mitocorrect_output/`                               |
| `-l`/`--logfile`         | The name of a file to which log messages will be written, saved to the output directory                                                                                                                                   | `mitocorrect.log`                                   |
| `-r`/`--detailedresults` | If a file name is supplied, mitocorrect will write a detailed table of individual potential annotation statistics and scores to filtering_results.tsv in the output directory                                             | *not written*                                       |
| `-1`/`--onefile`         | If a file name is supplied, mitocorrect will write all corrected sequences to a single genbank-formate file with this name in the output directory                                                                        | *each input file written to a separate output file* |

#### Uncommon optional arguments

| Argument                   | Description                                                                                                                                                                                                                 | Default                                                                                                          |
|----------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------|
| `-n`/`--namevariants`      | The path to a file specifying variant gene names to dd to those retrieved from the list [here](https://github.com/tjcreedy/genenames), in the same format                                                                   | *none*                                                                                                           |
| `-m`/`--maxinternalstops`  | The number of stops permitted in the translation of annotations                                                                                                                                                             | `0`                                                                                                              |
| `-e`/`--framefree`         | If used, mitocorrect wil allow codon search strings to not start in frame with one another; *not generally recommended*                                                                                                     | *N/A*                                                                                                            |
| `-p`/`--potentialfeatures` | If used, mitocorrect will add all potential annotations as features to the output genbank-format file(s); *not generally recommended*                                                                                       | *only selected annotations recorded*                                                                             |
| `-f`/`--fullassessment`    | If used, mitocorrect will assess the alignment of all potential annotations against the profile, even if their overlap score means they cannot possibly be selected; *not generally recommend, causes substantial slowdown* | *only align and score potential annotations that have good enough overlap scores to be contenders for selection* |

### Details

### Scoring and selection
Mitocorrect selects a potential annotation (specifically, a start/stop codon combination) based 
on, effectively, a set of nested thresholds:
1. Any start/stop codon combinations not in frame with one another are rejected (unless
`-e`/`--framefree` is set)
2. Any start/stop codon combinations where the translation of the potential annotation region, 
excluding any final stop, has more stops than specified by `-m`/`--maxinternalstops` are 
rejected
3. Mitocorrect calculates positionscore, the sum of the absolute distances in base pairs of each
of the start/stop codons for a given potential annotation from the overlap position given in the 
specifications, multiplied by positionweight for that end and gene
4. The potential annotation with the lowest (best) score is aligned against the profile and two
further scores are computed:
   1. the alignscore, the sum of the absolute distance of the aligned sequence start/stop from
the modal start/stop of the profile sequences, multiplied by the alignweight
   2. the indelscore, the number of nucleotide bases representing insertions or deletions in the
first and second half of the sequence of the potential annotation relative to the consensus are 
counted and multiplied by the indelweight 
   3. These three scores are summed to form the final score for the potential annotation 
5. Any subsequent potential annotation with a positionscore greater than the lowest current total
score is rejected without performing alignment, as negative scores are not possible (unless 
`-f`/`-fullassessment` is set). Any that are not rejected are aligned and recieve an alignscore,
indelscore and total score.

The potential annotation with the lowest total score is selected.

### To do

* Write a script that assesses a set of mitogenomes and outputs a specification table.