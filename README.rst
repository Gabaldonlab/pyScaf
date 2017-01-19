.. contents:: Table of Contents

pyScaf
======

pyScaf orders contigs from genome assemblies utilising several types of information:

- paired-end (PE) and/or mate-pair libraries ([NGS-based mode](#ngs-based-scaffolding))
- long reads ([NGS-based mode](#scaffolding-based-on-long-reads))
- synteny to the genome of some related species ([reference-based mode](#reference-based-scaffolding))

=================
Scaffolding modes
=================

NGS-based scaffolding
~~~~~~~~~~~~~~~~~~~~~
This is under development... Stay tuned. 

Scaffolding based on long reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Experimental version available.

Reference-based scaffolding
~~~~~~~~~~~~~~~~~~~~~~~~~~~
In reference-based mode, pyScaf uses synteny to the genome of closely related species in order to order contigs and estimate distances between adjacent contigs.

Contigs are aligned globally (end-to-end) onto reference chromosomes, ignoring:

- matches not satisfying cut-offs (`--identity` and `--overlap`)
- suboptimal matches (only best match of each query to reference is kept) 
- and removing overlapping matches on reference. 

In preliminary tests, pyScaf performed superbly on simulated heterozygous genomes based on *C. parapsilosis* (13 Mb; CANPA) and *A. thaliana* (119 Mb; ARATH) chromosomes, reconstructing correctly all chromosomes always for CANPA and nearly always for ARATH (`Figures in dropbox <https://www.dropbox.com/sh/bb7lwggo40xrwtc/AAAZ7pByVQQQ-WhUXZVeJaZVa/pyScaf?dl=0>`_, `CANPA table <https://docs.google.com/spreadsheets/d/1InBExy-qKDLj-upd8tlPItVSKc4mLepZjZxB31ii9OY/edit#gid=2036953672>`_, `ARATH table <https://docs.google.com/spreadsheets/d/1InBExy-qKDLj-upd8tlPItVSKc4mLepZjZxB31ii9OY/edit#gid=1920757821>`_).  
Runs took ~0.5 min for CANPA on `4 CPUs` and ~2 min for ARATH on `16 CPUs`. 

**Important remarks:**

- Reduce your assembly before (fasta2homozygous.py) as any redundancy will likely break the synteny.
- pyScaf works better with contigs than scaffolds, as scaffolds are often affected by mis-assemblies (no *de novo assembler* / scaffolder is perfect...), which breaks synteny. 
- pyScaf works very well if divergence between reference genome and assembled contigs is below 20% at nucleotide level. 
- pyScaf deals with large rearrangements ie. deletions, insertion, inversions, translocations. **Note however, this is experimental implementation!**
- Consider closing gaps after scaffolding. 

=====
Usage
=====
Dependencies
~~~~~~~~~~~~
- `LAST v700+ <http://last.cbrc.jp/>`_
- `FastaIndex <https://github.com/lpryszcz/FastaIndex>`_

Parameters
~~~~~~~~~~
Given reference genome, the program generates pairwise genome alignment (dotplots) by default. 

- Genral options:

  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        assembly FASTA file
  -o OUTPUT, --output OUTPUT
                        output stream [scaffolds.fa]
  -t THREADS, --threads THREADS
                        max no. of threads to run [4]
  --log LOG             output log to [stderr]
  --dotplot
                        generate dotplot as [png]
  --version             show program's version number and exit

- Reference-based scaffolding options:

  -r REF, --ref REF, --reference REF
                        reference FastA file
  --identity IDENTITY   min. identity [0.33]
  --overlap OVERLAP     min. overlap  [0.66]
  -g MAXGAP, --maxgap MAXGAP
                        max. distance between adjacent contigs [0.01 * assembly_size]
  --norearrangements    high identity mode (rearrangements not allowed)

- Long read-based scaffolding options (EXPERIMENTAL!): 

  -n LONGREADS, --longreads LONGREADS
                        FastQ/FastA file(s) with PacBio/ONT reads

- NGS-based scaffolding options (!NOT IMPLEMENTED!):

  -i FASTQ, --fastq FASTQ
                        FASTQ PE/MP files
  -j JOINS, --joins JOINS
                        min pairs to join contigs [5]
  -a LINKRATIO, --linkratio LINKRATIO
                        max link ratio between two best contig pairs [0.7]
  -l LOAD, --load LOAD  align subset of reads [0.2]
  -q MAPQ, --mapq MAPQ  min mapping quality [10]


Test run
~~~~~~~~
To perform reference-based assembly, provide assembled contigs and reference genome in FastA format.
Dotplots of below runs can be found in [docs](/docs).
If you wish to skip dotplot generation (ie. no X11 on your system), provide `--dotplot ''` parameter. 

.. code-block:: bash

    # scaffold homogenised assembly (reduced contigs)
    ./pyScaf.py -f test/contigs.reduced.fa -r test/ref.fa -o test/contigs.reduced.ref.fa

    # scaffold reduced contigs using global mode (no norearrangements allowed)
    ./pyScaf.py -f test/contigs.reduced.fa -r test/ref.fa -o test/contigs.reduced.ref.global.fa --norearrangements

    # scaffold heterozygous assembly (de novo assembled contigs)
    ./pyScaf.py -f test/contigs.fa -r test/ref.fa -o test/contigs.ref.fa

    # scaffold reduced contigs using long reads
    ## pacbio
    ./pyScaf.py -f test/contigs.reduced.fa -n test/pacbio.fq.gz -o test/contigs.reduced.pacbio.fa
    ## nanopore
    ./pyScaf.py -f test/contigs.reduced.fa -n test/nanopore.fa.gz -o test/contigs.reduced.nanopore.fa

    # generate dotplot
    lastdb test/ref.fa
    lastal -f TAB test/ref.fa test/contigs.reduced.pacbio.fa | last-dotplot - test/contigs.reduced.pacbio.fa.ref.png
    lastal -f TAB test/ref.fa test/contigs.reduced.nanopore.fa | last-dotplot - test/contigs.reduced.nanopore.fa.ref.png

    # clean-up
    #rm test/contigs.{,reduced.}fa.* test/ref.fa.* test/*.{nanopore,pacbio,ref}* test/*.log


================
Proof of concept
================
pyScaf is under heavy development right now.
Nevertheless, the reference-based mode is functional and produces meaningful assemblies. Moverover, it has been implemented in `Redundans <https://github.com/lpryszcz/redundans>`_.

For more info, have a look in `workbook <https://docs.google.com/document/d/1WNw6FYZXNI2sKJ1hBZ0LI9CWJSQ-BTQID7jL9lLvYaA/edit?usp=sharing>`_. 

