# External cohorts for species-level replication (class X): search of 2026-11-11

**What qualifies.** Faecal **shotgun** metagenomes and a faecal **metabolome** from the
**same stool samples**, CRC and controls, and data we can obtain. Ideally the metabolomics
platform matches Yachida's (CE-TOFMS).

| Study | Paired shotgun + metabolome | n, groups | Platform | Data | Verdict |
|---|---|---|---|---|---|
| **Nofinska et al., *Microbiome* 2026** ([PMC13732858](https://pmc.ncbi.nlm.nih.gov/articles/PMC13732858)) | **Yes**, a subset of 270 (102 early-onset CRC, 54 young controls, 68 late-onset CRC, 46 older controls); stages 0 to IV | 510 enrolled, 324 matched | **CE-TOFMS** (same as Yachida) + HiSeq2500 shotgun | Metagenomes: DDBJ **PRJDB19991, PRJDB37508**. Metabolome: **not deposited**; the data statement does not mention it. | **Best candidate, with two blockers.** (1) It comes from the same NCC group as Yachida 2019, cites the same collection protocol, and "S.Y. provided the raw data", so **participant overlap with Yachida is likely** and has to be ruled out sample by sample against the BioSample accessions. (2) The metabolite table must be requested from the co-corresponding authors (S. Yachida, U. Osaka; T. Yamada, Inst. Science Tokyo). |
| Sun et al., *Cancer Cell* 2024 (doi 10.1016/j.ccell.2024.07.005) | Unknown: the publisher pages blocked automated reading | 1,251 across 4 cohorts (CRC, adenoma, controls), plasma and faecal metabolomics | Unknown | Search snippets mention CNGB CNP0005876 / CNP0005882; **unverified** | Check by hand: Data and code availability section. |
| Borenstein curated collection v2.1 ([wiki](https://github.com/borenstein-lab/microbiome-metabolome-curated-data/wiki/)) | 14 paired datasets | | | Public | The only shotgun CRC dataset is **Yachida**. Sinha and Kim are 16S; Erawijantari is already used. **Nothing new.** |
| Spathakis et al., *Pharmaceuticals* 2026 ([PMC13708536](https://pmc.ncbi.nlm.nih.gov/articles/PMC13708536)) | **No**, by its own statement: "metabolomics datasets were not paired with the metagenomic samples" | 1,534 metagenomes (curatedMetagenomicData) + Metabolomics Workbench ST002787 (30/37/35) | | Public | Unpaired. It cannot replicate a species–metabolite *association*, only a metabolite direction. |
| Zhang et al., *Cancer Med* 2023 ([PMC10417192](https://pmc.ncbi.nlm.nih.gov/articles/PMC10417192)) | No: shotgun on 14 samples, no metabolomics | 50 CRC / 52 controls (16S) | | On request | Not usable. |

## Next steps (none taken)
1. **Nofinska overlap check. Accessions: done, 2026-11-11.** Through ENA, PRJDB19991 has 40
   samples and PRJDB37508 has 470 (510, the number enrolled). **None of them shares a
   BioSample accession with Yachida's 347.** That rules out the same deposited samples. It
   does **not** rule out the same people re-sampled or re-deposited: same group, same
   ethics approval (2013-244), same protocol. Subject-level confirmation needs their
   Table S1 or the authors.
2. If enough independent samples remain, **draft a data request** to the corresponding
   authors for the CE-TOFMS table (to be sent only after the user approves the text).
   The metagenomes would also need profiling against the same GTDB reference as Yachida.
   That is a large download, so state its GB and get a yes first.
3. Sun et al. 2024: read the Data and code availability section by hand.
