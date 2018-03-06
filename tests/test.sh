Rscript blatCouple.R tests/data/testSeq-1.R2.psl.gz tests/data/testSeq-1.R1.psl.gz \
    -k tests/data/testSeq-1.R2.key.csv tests/data/testSeq-1.R1.key.csv \
    -o tests/testSeq-1.uniq.csv --condSites tests/testSeq-1.cond.csv \
    --chimera tests/testSeq-1.chimera.rds --multihit tests/testSeq-1.multihit.rds \
    --refGenome hg38
