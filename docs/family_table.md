# Family table for clinicians

PGA needs a pedigree (PED), but nobody should have to write one by hand. The
clinician fills in a simple table with words, PGA draws the family and the
analysis only runs after the clinician approves the drawing.

```bash
pga family template --out family.tsv --lang en     # or --lang es
# fill it in with Excel / LibreOffice / any editor
pga family check family.tsv                         # drawing + table as understood + PED
pga family approve family.tsv --by "Dr Name"        # writes family.ped + family.approval.json
pga run-family --family family.tsv --reference GRCh38.fa --region chr2:... [--hpo ...]
```

## The table

One row per relative, sequenced or not. Relations are written relative to the
proband: proband, father, mother, brother, sister, son, daughter,
paternal/maternal grandfather or grandmother, paternal/maternal uncle or aunt,
or `other` (then give the ids of the parents in the `father` and `mother`
columns). Sex: male, female, unknown. Status: affected, unaffected, unknown.
`bam`: the file of a sequenced individual, relative to the table. Spanish words
are accepted (probando, padre, abuela materna, afectado, sano…).

One header line states how the proband's parents are related: unrelated,
first cousins, second cousins, uncle-niece, unknown. PGA adds the
grandparents and shared ancestors needed to close the loop, marks them
"(added)", and checks that the expected inbreeding coefficient of the
proband matches (1/16 for first cousins, 1/64 for second cousins, 1/8 for
uncle–niece). For first cousins it assumes the paternal grandfather and the
maternal grandmother are siblings; autosomal kinship does not depend on this
choice, but for X-linked analysis give the exact links with `father`/`mother`.

## Review and approval

`pga family check` writes `family.check.html`: the pedigree drawn with
[Pedixplorer](https://bioconductor.org/packages/Pedixplorer) (Bioconductor,
Artistic-2.0, successor of kinship2; R environment `pga env bootstrap pga-r`),
the table as PGA understood it, and the PED that will be used.

`pga family approve` is refused unless that exact version of the table has
been reviewed. It writes the PED and `family.approval.json` (SHA-256 of the
table and the PED, who approved and when). `pga run-family --family` refuses
to run if the table is not approved or has changed since; the approval is
recorded in the session.
