# PublicGenomicAgent — verified facts

Every figure below was measured on the date given, with the command next
to it, so anyone can reproduce it. Re-run the commands on the day this
sheet is shared and update the dates and numbers.

## What it does

Local workflow for families with rare Mendelian disease, aimed at
consanguineous populations. Input: confirmed PED, one BAM per sequenced
individual (panel or exome), optional HPO terms. Order: per-sample QC
gate → family verification → joint calling → rare-variant filter →
segregation (recessive homozygous, compound heterozygous with phase) →
prioritisation → self-contained HTML report. Research use only.

## Measured results

| Claim | Value | Date | How to reproduce |
|---|---|---|---|
| Unit tests in a clean environment (as GitHub CI) | 257 passed, 5 skipped | 2026-11-10 | `pytest tests/unit -q` |
| Unit + integration tests on the author's laptop | 262 unit + 43 integration passed (before E3) | 2026-11-08 | `pytest tests -q` |
| Family verification, true trio (GIAB HG002/3/4 exomes, chr21–22) | never fails; 98.6 % passed, 1.4 % warning at 200 sites; 100 % passed at ≥ 500 | 2026-10-31 | `scripts/calibrate_relatedness_giab.py`; [`relatedness_calibration.md`](relatedness_calibration.md) |
| Family verification, false father (HG005 declared father of HG002) | never passes; 100 % failed at ≥ 150 sites | 2026-10-31 | same |
| Below ~200 informative sites | "not verifiable", never "passed" | 2026-10-31 | same |
| E2 synthetic family on real GRCh38 (CYP2D6) | 8,743 gnomAD SNVs with matching REF, 0 discordant | 2026-11-08 | `scripts/make_e2_case.py --workdir DIR` |
| E2 rare-variant filter (`--max-af 0.01`) | 28 common variants removed | 2026-11-08 | same |
| E2 result | 1 recessive candidate, the seeded one (chr22:42129004 C>A, gnomAD AF 6.9e-7), ranked 1 | 2026-11-08 | same |
| Report opens without network | 0 external requests, 0 JavaScript errors in Chromium with all network blocked | 2026-11-08 | `pytest tests/integration/test_report_viewer.py` |

## Network and data

- The only step that uses the network is the gnomAD frequency lookup
  (`--max-af`). It sends region coordinates, not genotypes or reads, to
  the public gnomAD API and caches the answers locally.
- A report generated with reads contains individual genomic data.
  `--no-reads` leaves them out; `--redact-paths` replaces local file
  paths with file names.
- The embedded viewer is JBrowse 2 (Apache-2.0), compiled once and stored
  in the repository with its SHA-256.

## What has not been shown

- No real patient family has been analysed.
- Family verification was calibrated on chromosomes 21–22 of GRCh37
  exomes with a 1000 Genomes site panel. No false mother was tested, and
  the power against a non-parent of the same population is calculated,
  not measured. A panel with Middle Eastern frequencies has not been built.
- Tested on gene panels and single regions. Exome-wide runs with many
  regions and autozygosity mapping are not done yet.
- The local-pangenome benchmark on a synthetic segmental duplication
  (`benchmark_duplications.md`) cannot be reproduced from the repository
  yet: its generator script and demo module are not versioned. It is not
  cited as a result until it is rebuilt with a versioned script, a test,
  and a graph built without the target variant.
- E2 is a simulated family on real reference sequence and real
  population frequencies. It shows that the steps work together; it is
  not a clinical validation.
