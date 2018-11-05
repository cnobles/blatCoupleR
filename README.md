# blatCoupleR
[![Build Status](https://travis-ci.org/cnobles/blatCoupleR.svg?branch=master)](https://travis-ci.org/cnobles/blatCoupleR)

Join paired-end sequence alignments from BLAST-Like Alignment Tool (blat) output (psl format).

## arguments
### positional arguments:
  **[anchorPSL]** Alignment file (psl format) to couple, component of paired-end sequence that is anchored to a biological phenomenon.
  
  **[adriftPSL]** Alignment file (psl format) to couple, other component of paired-end sequnece that is generated by non-natural forces (shearing).

### optional arguments:
  **[-h, --help]** show help message and exit
  
  **[-k, --keys]** Key files which denote the consolidation of reads to unique sequences. See README for requirements. Mulitple key files should be input in the same order as their respective psl files.
  
  **[-o, --uniqOutput]** Output file for unique alignments. File types supported: .rds, .RData, .csv, and .tsv
  
  **[--condSites]** Output file for condensed sites, based on anchor alignments. Read counts and unique alignment length counts are reported for each unique site. Same file types supported as uniqOutput.
  
  **[--chimeras]** Output file for chimeric alignments. Same file types supported as uniqOutput.
  
  **[--multihits]** Output file for multihit alignments. Same file types supported as uniqOutput.
  
  **[--stat]** File name to be written in output directory of read couts for each sample. CSV file format. ie. test.stat.csv.

  **[-g, --refGenome]** Reference genome, needs to be installed through BSgenome (BioConductor).
  
  **[--maxAlignStart]** Maximum allowable distance from the start of the sequence to keep the alignment. Default = 5.
  
  **[--minPercentIdentity]** Minimal global (whole sequence) percent identity required to keep alignment. Default = 95 (0-100).
  
  **[--minTempLength]** Minimum value for paired template length to consider. Default = 30 (bps).
  
  **[--maxTempLength]** Maximum value for paired template length to consider. Default = 2500 (bps).
  
  **[--keepAltChr]** By default, blatCoupleR will remove alignments from psl files aligning to alternative chromosome sequences, ex. chr7_*_alt. Using this option will keep these alignments, which may increase multihit outputs.
  
  **[--readNamePattern]** Regular expression for pattern matching read names. Should not contain R1/R2/I1/I2 specific components. Default is [\w:-]+
  
  **[--saveImage]** Output file name for saved image. Include '.RData' (ie. debug.RData).

## dependencies
blatCouplR is coded in R and was developed on v3.4.0, though it should run with earlier versions given the appropriate dependencies. The script uses 9 additional packages:
  * argparse
  * pander
  * yaml
  * stringr
  * GenomicRanges
  * igraph
  * data.table
  * Matrix
  * BSgenome
  * BSgenome.Hsapiens.UCSC.hg38 (for tests)
