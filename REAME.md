### Table of Contents
- **[Scaffolding modes](#scaffolding-modes)**
  - **[NGS-based scaffolding](#ngs-based-scaffolding)**
  - **[Scaffolding based on long reads](#scaffolding-based-on-long-reads)
  - **[Reference-based scaffolding](#reference-based-scaffolding)**
- **[Usage](#usage)**

# pyScaf

pyScaf orders contigs from genome assemblies utilising several types of information:
- paired-end (PE) and/or mate-pair libraries ([NGS-based mode](#ngs-based-scaffolding))
- long reads ([NGS-based mode](#scaffolding-based-on-long-reads))
- synteny to the genome of some related species ([reference-based mode](#reference-based-scaffolding))

## Scaffolding modes

### NGS-based scaffolding
This is under development... Stay tuned. 

### Scaffolding based on long reads
This is under development... Stay tuned.

### Reference-based scaffolding
In reference-based mode, pyScaf uses synteny to the genome of closely related species in order to order contigs and estimate distances between adjacent contigs.

Contigs are aligned globally (end-to-end) onto reference chromosomes, ignoring:
- matches not satisfying cut-offs (`--identity` and `--overlap`)
- suboptimal matches (only best match of each query to reference is kept) 
- and removing overlapping matches on reference. 

In preliminary tests, pyScaf performed superbly on simulated heterozygous genomes based on *C. parapsilosis* (13 Mb; CANPA) and *A. thaliana* (119 Mb; ARATH) chromosomes, reconstructing correctly all chromosomes always for CANPA and nearly always for ARATH ([Figures in dropbox](https://www.dropbox.com/sh/bb7lwggo40xrwtc/AAAZ7pByVQQQ-WhUXZVeJaZVa/pyScaf?dl=0), [CANPA table](https://docs.google.com/spreadsheets/d/1InBExy-qKDLj-upd8tlPItVSKc4mLepZjZxB31ii9OY/edit#gid=2036953672), [ARATH table](https://docs.google.com/spreadsheets/d/1InBExy-qKDLj-upd8tlPItVSKc4mLepZjZxB31ii9OY/edit#gid=1920757821)).  
Runs took ~0.5 min for CANPA on `4 CPUs` and ~2 min for ARATH on `16 CPUs`. 

**Important remarks:**
- Reduce your assembly before (fasta2homozygous.py) as any redundancy will likely break the synteny.
- pyScaf works better with contigs than scaffolds, as scaffolds are often affected by mis-assemblies (no *de novo assembler* / scaffolder is perfect...), which breaks synteny. 
- pyScaf works very well if divergence between reference genome and assembled contigs is below 20% at nucleotide level. 
- pyScaf deals with large rearrangements ie. deletions, insertion, inversions, translocations. **Note however, this is experimental implementation!**
- Consider closing gaps after scaffolding. 

# Usage

```bash
./pyScaf.py -f test/run1/contigs.reduced.fa -r test/ref.fa
```

## Proof-of-concept


