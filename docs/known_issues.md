# Bugs conocidos y limitaciones

Este documento recoge bugs conocidos, limitaciones actuales y
trabajo pendiente de PublicGenomicAgent.

## Bugs activos

### Simulador de variantes: solo aplica a lecturas que cubren
### completamente el REF de la variante

**Fichero**: `tests/fixtures/make_bam_with_variant.py`

**Síntoma**: con variantes cuyo REF es largo (por ejemplo, un
homopolímero de 20 T's), solo una fracción muy pequeña de las
lecturas (0.24% con cobertura 30x) porta la variante. `bcftools
call` no detecta la variante por falta de evidencia.

**Causa**: la función `build_read` exige que la lectura cubra
**completamente** el REF de la variante:

    covers_variant = (start_idx <= ref_idx) and
                     (ref_idx + var_len_ref <= read_end_ref)

Con REF de 20 bp y lecturas de 150 bp, la ventana válida de
`start` es muy estrecha.

**Workaround actual**: reducir el REF de la variante en el VCF a
una subsecuencia corta (por ejemplo, 5 T's de las 20) y ejecutar
el simulador con esa variante reducida. Funciona para validar el
flujo, pero pierde información sobre el tamaño real del indel.

**Fix pendiente**: modificar `build_read` para aplicar la variante
a cualquier lectura que **solape parcialmente** con la región del
REF. Para delecciones largas, la CIGAR puede tener M pequeños antes
o después del `D`.

### Simulador de variantes: no soporta múltiples variantes por BAM

**Fichero**: `tests/fixtures/make_bam_with_variant.py`

**Síntoma**: un individuo heterocigoto compuesto (dos alelos
distintos en el mismo locus, como el hijo KSA001 real) no se puede
modelar con un solo BAM.

**Workaround actual**: modelar el hijo como heterocigoto simple
(un alelo alt). Simplificación documentada en `demo_ksa001_trio.md`.

**Fix pendiente**: aceptar una lista de variantes (o un VCF
completo) y aplicar cada una a las lecturas que la cubran.

## Limitaciones conocidas

### Extracción HPO (RDMA) requiere GPU ≥20 GB

Ver `docs/phenotype_layer.md` sección 14.

### No alineamos

PublicGenomicAgent empieza con BAMs ya alineados a hg38. El
alineamiento (FASTQ → BAM) no es parte del sistema.

### LIRICAL requiere Java 17+

Y tiene un bug de duplicados en HGNC/NCBI que se arregla con el
fix documentado en `docs/lirical_setup.md`.

## Trabajo pendiente

- Arreglar los dos bugs del simulador de variantes.
- Soporte para descargar BAMs reales de SRA/ENA.
- Integración con cohortes acumulativas (`docs/local_pangenome.md`
  sección 14).
- `pdf_extract` para cargar PDFs clínicos.
- `annotate` con VEP/snpEff.
- Crash test de instalación desde clon limpio.
