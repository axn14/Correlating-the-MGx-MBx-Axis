# Microbial sources of CRC-associated faecal metabolites

This wiki documents version 2.0 of the analysis. Pages marked "v1.0-thesis" describe the
earlier pipeline and are kept for reference only.

The question is which gut bacteria produce the faecal metabolites that distinguish colorectal
cancer (CRC) from health. The data are public: 277 paired faecal metagenomes and metabolomes
from Yachida et al. (2019) for discovery, and 42 post-gastrectomy samples from Erawijantari
et al. (2020), on the same metabolomics platform, for replication.

## Pages

| Page | What it covers |
|---|---|
| [[What changed in v2.0]] | Why v2 replaces the thesis pipeline, and which earlier results no longer hold |
| [[Evidence classes]] | The five independence classes, the tier rule and the gates |
| [[Key results]] | The v2 numbers |
| [[Pipeline Overview]], [[Methods]], [[Evidence Streams]], [[Data Sources]], [[Installation]] | v1.0-thesis documentation, kept for reference |

## The answer in brief

Counting each class of evidence once, 82 of 4,293 species–metabolite pairs reach the
Supported tier or better. A permutation null expects about 45 by chance. For the 30
CRC-associated metabolites, Supported links do not exceed chance. The standard does recover
known producers where gene evidence points at a single species: *Ruminococcus gnavus* and
tryptamine, and *Veillonella atypica* and putrescine. The rules for all of this are in
`v2/CONTRACTS.md`. Each one was committed before the results it governs.
