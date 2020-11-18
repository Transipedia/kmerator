# kmerator

## A prototype for kmer decomposition and extraction of transcript or gene-specific k-mers

Kmerator is a prototype tool designed for prediction of specific k-mers (named tags) from an input sequence, considering the sequence of a reference genome and transcriptome. Kmerator use jellyfish to create two requestable indexes from these references, decompose your input sequence and count the occurences of each kmers of the input in the genome and transcriptome. The occurrences counts are then interpreted, and passed to different filters to select tags from the input.

Kmerator strictly depends of an Ensembl Annotation, you can find it there : https://www.ensembl.org/info/data/ftp/index.html

You’ll need a genome and a transcriptome, fasta format. You can reuse the jellyfish of the genome in the options. 
For a more complete k-mer filtering, we advice to merge the coding (cDNA) and the non-coding (ncRNA) as one unique transcript reference.

You can choose different k-mer filters with the “level” option : 
- gene 
- transcript 
- chimeras (must not be fount in your reference, genome or transcriptome)

If you want to use kmerator on unannotated RNA signatures, you have to use the -u option. In this case only "transcript" and "chimeras" levels are available.

## Dependencies :

- Kmerator is a julia script : Julia >= v1.4.2 (https://julialang.org/downloads/).
- kmerator ask to install missing packages automatically so you may not install ParallelDataTransfer, ArgParse, FastaIO, JSON and HTTP packages.
- Jellyfish >= v2.0


## Usage
```
kmerator.jl [--selection [SELECTION...]] [--appris APPRIS] [-u]
                   [-s] [-v] -g GENOME -t TRANSCRIPTOME -l LEVEL
                   [-o OUTPUT] [--length LENGTH] -p [nb-threads]
                   [--threshold THRESHOLD] fasta_file
```

## arguments
```
  --selection [SELECTION...]
                        list of gene to select directly inside your
                        fasta transcriptome file (ENST, ENSG, official_symbol)
  --appris APPRIS        indicate :
                        'homo_sapiens','mus_musculus','rattus_norvegicus','danio_rerio','sus_scrofa',''
                        or virtually any specie available on APPRIS
                        database (Rodriguez JM, Rodriguez-Rivas J,
                        Domenico TD, Vázquez J, Valencia A, and Tress
                        ML. Nucleic Acids Res. Database issue; 2017
                        Oct 23.).     A gene have multiple possible
                        sequence depending of its variants. this
                        option select the principal transcript of a
                        gene, based of APPRIS database, or by length
                        if no data or no internet connection. If 'APPRIS’ option
                        and not available data or no principal
                        transcript (non-coding genes), use length
                        instead
  -u, --unannotated     activated if the provided initial fasta file
                        correspond to an annotation external from
                        Ensembl. 
  -s, --stringent       FOR GENE LEVEL ONLY : If you want select tags present 
                        in ALL variants of the corresponding gene 
  -v, --verbose         If you want this script talk too much
  -g, --genome GENOME   the genome fasta (.fa) or index by jellyfish
                        for kmer request
  -t, --transcriptome TRANSCRIPTOME
                        The transcriptome fasta (.fa) (FASTA ONLY) for
                        kmer request and transcriptional variant
                        informations
  -l, --level LEVEL     Type 'gene', 'transcript' or 'chimera' to
                        extract tags at these different
                        levels. 'chimera' option must be done with
                        'unannotated' option!.
  -o, --output OUTPUT   directory of output (default: ".")
  --length LENGTH       length required for the kmer generation (type:
                        Int64, default: 31)
  --threshold THRESHOLD
                        FOR GENE LEVEL ONLY : Minimum fraction of
                        annotated transcripts containing this kmer to
                        admit a tag (default: 0.5)
  -p, --procs PROCS     Run n local processes
  -h, --help            show this help message and exit

```
