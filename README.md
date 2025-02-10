# Kmerator

## Prototype for decomposition of transcript or gene sequences and extraction of their specific k-mers

ref: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8221386/>

Kmerator is a prototype tool designed for the prediction of specific k-mers (also called tags) from input sequences, considering a reference genome and an ENSEMBL-like transcriptome. From these specific k-mers, it also outputs their corresponding specific contigs which are sequences of consecutive k-mers (overlapping length between k-mers must be k-1, otherwise, it's a new contig). You need to provide kmerator with a jellifsh index of the reference genome. Kmerator itself builds a jellyfish index of the reference transcriptome (by default the latest available version of Ensembl). It then  decomposes your input transcript or gene sequences to count the occurences of each k-mer in the genome and transcriptome. Number of occurrences are then interpreted, in different manners, to select specific k-mer from your input. 

Before using kmerator, a jellyfish index of the reference genome must be created. kmerator automatically creates a dataset according to the species and the desired release number (by default, homo_sapiens and the latest version). The dataset is composed of 4 files per species/version: a jellyfish index of the modified transcriptome (cDNA + ncRNA - alternative chormosomes) from Ensembl, a binary file representing the same transcriptome, another binary file containing general information on the genes of the transcriptome and a report file.


#### Specific kmers

![](https://github.com/Transipedia/kmerator/raw/main/img/specific-kmers.png)

#### Specific contigs

![](https://github.com/Transipedia/kmerator/raw/main/img/specific-contigs.png)

## Dependencies

- Python >= v3.7
- Jellyfish >= 2.0


## Installation

### Option 1 (preferred)

Install with pip

```
pip3 install kmerator
```


With pip, it is also easy to install an older version

```
pip install kmerator==1.0.0
```

### Option 2

Installation from github

```
git clone https://github.com/Transipedia/kmerator.git
ln -s $PWD/kmerator3/kmerator/kmerator.py /usr/local/bin/kmerator  # or somewhere in your $PATH
```


## How to use kmerator

Before all, remember that kmerator needs a jellyfish index of the genome. You must build it according to the species you are studying. You can store and name the index file whatever you want. Please note that you must **use the jellyfish -C option** when building the reference genome index.

### Configuration file

The arguments to run kmerator are numerous, so to reduce the number of arguments to enter, it is advisable to edit the configuration file with the command :

```
kmerator -e
```

By filling in the `datadir` and `genome` directives, you will avoid having to re-enter the `--datadir` and `--genome` arguments systematically. If you are working on a species other than Human, you can also fill in the specie directive. And in a long-term project, you may want to set a release number.

### Execute requests

There are two main cases:

- you find for specific k-mers for annotated genes or transcripts : use the `--selection` option, followed by:
	- the list of gene and/or transcripts separated by a space
	- or a file with the list of genes/transcripts. Separator could by a space, a tab or a newline, and comments are allowed (`#`)
- you find for specific k-mers of unannotated sequences : use the `--fasta-file` option, followed by a fasta file containing yours requests. In case of you focuses on chimeras, add the `--chimera` option

**Examples:**

```
kmerator -s npm1 brca2 ENST00000255409 ENSG00000159216    # you can mix genes and transcripts
kmerator -s genes.txt                                     # you can also use a file with gene list
kmerator -f file.fa                                       # give a fasta file fr unannotated sequences
```

**Note** the above commands assume that the configuration file contains at least the `datadir` and `genome` directives, the default species is homo_sapiens and the last available version will be used (if it is not present in datadir, kmerator will propose the construction of a dataset automatically)


### Note the difference between genes and transcripts

- When you are looking for specific kmers of a **gene** (symbol, alias or Ensembl name), kmerator fetch sequence of its canonical transcript, extracts kmers and keep those that found only in the gene.
- When you are looking for a **transcript**, kmerator only keeps the kmer found in the transcript, and only in that transcript. If isoforms completely cover the transcript, no kmer will be kept.

### Datasets

To work, kmerator needs a jellyfish index of the genome, a jellyfish index of the transcriptome and various files. You will have to make the jellyfish genome index manually. Instead, kmerator builds the jellyfish transcriptome index and the files it needs, which we call datasets. There is one dataset per species and per transcriptome version. When kemrator does not find (in datadir) the requested transcriptome release (by default, the latest available on Ensembl), it offers to automatically build the dataset in question. In addition, dataset management options are available:

```
kmerator -l            # list local datasets
kmerator -u            # find last release on Ensembl, and build dataset if not present
kmerator --mk-dataset  # build dataset according to -r <release> and -S <specie> arguments
kmerator --rm-dataset  # delete dataset according to -r <release> and -S <specie> arguments
```

### Info

You can get information about gene, using the `--info` parameters. Like previous argument, you can enter a mix of symbol gene, Ensembl gene or Ensembl transcript (ex: braf, ENSG00000157764, ENST00000646427), or even use a text file with the list of requested genes separated by space, tab or newline (comments: '#').
In addition, you can add `--all` argument to get extended information, like transcript sequences.

```
kmerator --info MLLT3 ENSG00000157764 # info about the MLLT3 and BRAF genes
kmerator --info mllt3 -r 109 -S mouse # info about mllt3 against mus_musculus release 109
kmerator --info MMLT3 --all           # extended info about the MMLT3 gene
kmerator --info genes.txt             # genes/transcripts are in a file
```

## All arguments

```
options:
  -h, --help            show this help message and exit
  -s SELECTION [SELECTION ...], --selection SELECTION [SELECTION ...] 
                        list of gene IDs (ENSG, gene Symbol or alias) or transcript 
                        IDs (ENST) from which you want to extract specific kmers from. 
                        For genes, kmerator search specific kmers along the gene. 
                        For transcripts, it search specific kmers to the transcript. 
                        You can also give a file with yours genes/transcripts 
                        separated by space, tab or newline. If you want to use 
                        your own unannotated sequences, you must give your fasta 
                        file with --fasta-file option.
  -f FASTA_FILE, --fasta-file FASTA_FILE
                        Use this option when yours sequences are unannonated or 
                        provided by a annotation file external from Ensembl. 
                        Otherwise, use --selection option.
  -d DATADIR, --datadir DATADIR
                        Storage directory for kmerator datasets.We recommend to 
                        set this parameter by editing the configuration file 
                        (kmerator --edit)
  -g GENOME, --genome GENOME
                        Genome jellyfish index (.jf) to use for k-mers requests.
  -S SPECIE, --specie SPECIE
                        indicate a specie referenced in Ensembl, to help, follow
                        the link https://rest.ensembl.org/documentation/info/species.
                        You can use the 'name', the 'display_name' or any 'alias'.
                        For example human or homo_sapiens are valid
                        (default:human).
  -k KMER_LENGTH, --kmer-length KMER_LENGTH               
                        k-mer length that you want to use (default 31).
  -r RELEASE, --release RELEASE         
                        release of transcriptome (default: last).
  --stringent           Only for genes with '-s/--selection' option: use this
                        option if you want to select gene-specific k-mers present
                        in ALL transcripts for your gene. If false, a k-mer is 
                        considered as gene-specific if present in at least one
                        isoform of your gene of interest.
  -T MAX_ON_TRANSCRIPTOME, --max-on-transcriptome MAX_ON_TRANSCRIPTOME
                        Only for genes with '-f/--fasta-file' option: with
                        unanotated data, specific kmers are not supposed to be
                        found in the transcriptome, you can change this behavior
                        for special cases (default: 0)
  -G MAX_ON_GENOME, --max-on-genome MAX_ON_GENOME        
                        Only for genes with '-f/--fasta-file' option: typically, 
                        specific kmers are not supposed to be found more than 
                        once in the genome, you can change
                        this behavior for special cases, like chimera (default: 1)
  -o OUTPUT, --output OUTPUT 
                        output directory, created if not exists (default: 'output')
  -t THREAD, --thread THREAD
                        run n process simultaneously (default: 1) 
  --tmpdir TMPDIR       directory to temporary file (default: /tmp/kmerator_<random>) 
  -D, --debug           Show more details while Kmerator is running.
  --keep                keep intermediate files (sequences, indexes, separate 
                        kmers and contigs files).
  -y, --yes             assumes 'yes' as the prompt answer, run non-interactively.
  -e, --edit-config     Edit config file 
  -l, --list-dataset, --list-datasets 
                        list the local datasets (based on the datadir option).
  --rm-dataset          remove a dataset, according with --specie and 
                        --release options
  --mk-dataset          make a dataset, according with --specie and --release options
  --last-avail, --last-available    
                        last release available on Ensembl
  -u, --update-dataset  builds a new dataset if a new version is found on Ensembl
  --info gene [gene ...] 
                        get some information about gene. Multiple entries are 
                        allowed or a text file with list of genes
  -a, --all             only with '--info' option. Give more info, like 
                        transcript sequences                       
  -v, --version         show program's version number and exit
```

## Output

kmerator generate 3 files:

* kmers.fa
* contigs.fa
* report.md

extract of kmers.fa:
```
>RUNX1:ENST00000675419.kmer3394 ct:3 tr:3/17
TGAAGAGTATTTGAAAGCAGGACTTCAGAAC
```
* `kmer3394`: the first base is at this position in the canonical transcript or the given sequence (1 based)
* `ct:3`: the kmer is associated to the 3rd contig of this gene/transcript, relative to `contigs.fa` file
* `ex3/17`: the kmer is located in the 3rd exon of the 17

extract of contigs.fa:
```
>RUNX1:ENST00000675419.contig_3 (at position 2314)
ACTTCTTTGGGCCTCATAAACAACCACAGAACCACAAGTTGGGTAGCCTGGCAGTGTCAGAAGTCTGAACCCAG
```
* `contig_3`: contig count, relative to `ct:3` of  `kmers.fa` headers
* `position 2314`: the first base of the contig is at this position (1 based)


## References

[1] Guillaume Marçais, Carl Kingsford, A fast, lock-free approach for efficient parallel counting of occurrences of k-mers, Bioinformatics, Volume 27, Issue 6, 15 March 2011, Pages 764–770, https://doi.org/10.1093/bioinformatics/btr011
[2] Rodriguez JM, et al. Nucleic Acids Res. Database issue; 2017 Oct 23
