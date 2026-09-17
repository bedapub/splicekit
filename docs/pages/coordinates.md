# Genomic coordinates

All genomic coordinates splicekit operates with are **0-based, left+right inclusive**. E.g. the range 100-103 includes coordinates 100, 101, 102 and 103; the first coordinate is 0.

More specifically:

- **Feature-specific**: all feature coordinates (junctions, anchors, exons) are given in numeric sort order regardless of strand — `feature_start` is always `<` `feature_stop`. Example: `chr1+_100_102` represents a junction spanning coordinates `[100, 101, 102]`; `chr1-_100_102` represents a junction spanning coordinates `[102, 101, 100]`.
- **Junction-specific**: junction coordinates cover/overlap 1 nucleotide of each adjoining exon.
    - `chr1+_100_200` — a junction on chromosome 1 (`+` strand) from `[100..200]`. 100 is the last nucleotide of the upstream exon, 200 is the first nucleotide of the downstream exon.
    - `chr1-_100_200` — a junction on chromosome 1 (`-` strand) from `[100..200]`. 200 is the last nucleotide of the upstream exon, 100 is the first nucleotide of the downstream exon.

!!! important
    RefSeq and Ensembl GTF files are 1-indexed. When splicekit reads files from RefSeq/Ensembl, it performs `coordinate -= 1` on every coordinate to keep them consistent with splicekit's internal 0-indexed structures.
