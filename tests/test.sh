#!/usr/bin/env bash
set -ev

# Test script for functionality
Rscript blatCouple.R tests/data/testSeq-1.R2.psl.gz tests/data/testSeq-1.R1.psl.gz \
    -k tests/data/testSeq-1.R2.key.csv tests/data/testSeq-1.R1.key.csv \
    -o tests/testSeq-1.uniq.csv --condSites tests/testSeq-1.cond.csv \
    --chimera tests/testSeq-1.chimera.rds --multihit tests/testSeq-1.multihit.rds \
    --refGenome hg38 --stat tests/testSeq-1.stat.csv

# Check output for correct findings, 50 reads / alignments and 5 sites
test_uniq_len=$(cat tests/testSeq-1.uniq.csv | sed '/seqnames/d' | wc -l)
test_cond_len=$(cat tests/testSeq-1.cond.csv | sed '/seqnames/d' | wc -l)

if [ ! $test_uniq_len = 50 ] | [ ! $test_cond_len = 5 ]; then
    exit 1
fi

head tests/testSeq-1.uniq.csv

head tests/testSeq-1.cond.csv

cat tests/testSeq-1.stat.csv

# Cleanup test directory
rm tests/testSeq-1.*

echo "Passed all tests."
#exit
