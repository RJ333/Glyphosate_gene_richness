#!/bin/sh

awk 'BEGIN { OFS = "\t" } ; {split ($1, a, "_"); split (a[3], b, "\\."); print $2, a[2], b[1]}' test_file.tsv