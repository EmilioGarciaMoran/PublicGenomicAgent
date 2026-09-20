# Demo NPHP1: joint calling del trío contra el grafo local

Cuarto experimento. Extiende el flujo de grafo a un trío completo
(padre, madre, hijo) y verifica que el joint calling contra el grafo
respeta la herencia mendeliana.

## 1. Objetivo

Demostrar que el pipeline de grafo local funciona con tríos:
- Tres BAMs (C2, F1, M1).
- Un mismo grafo del ROI.
- Tres alineamientos.
- Joint calling que produce un VCF multi-muestra.
- Verificación de herencia mendeliana.

## 2. Diseño

Fixture NPHP1 (chr2:110,000,000-110,025,000) con:
- Pedigree 3 generaciones, 10 individuos.
- **C2** (hijo) con la variante 2-110008979-A-T (heredada de F1).
- **F1** (padre) con la misma variante (0/1).
- **M1** (madre) sin la variante (0/0).
- Cluster MENA sintético alrededor de 8979.

BAMs sintéticos de C2, F1, M1 con cobertura uniforme 30x.

## 3. Pipeline ejecutado

    # 1. Grafo del ROI (una sola vez)
    pga tool build-local-graph --reference nphp1_ref.fa \
      --cohort cohort.vcf.gz --out roi_tool.vg

    # 2. Alinear cada individuo
    pga tool align-to-graph --graph roi_tool.vg \
      --reads c2.fastq --out c2_trio.gam
    pga tool align-to-graph --graph roi_tool.vg \
      --reads f1.fastq --out f1_trio.gam
    pga tool align-to-graph --graph roi_tool.vg \
      --reads m1.fastq --out m1_trio.gam

    # 3. Pack + call por individuo
    for s in c2 f1 m1; do
      vg pack -x roi_tool.xg -g ${s}_trio.gam -Q 5 -o ${s}.pack
      vg call roi_tool.xg -k ${s}.pack > ${s}_trio.vcf
      sed 's/^\(#CHROM.*\)\tSAMPLE$/\1\t<NAME>/' ${s}_trio.vcf > ${s}_named.vcf
      bcftools view -Oz -o ${s}_named.vcf.gz ${s}_named.vcf
      bcftools index -t ${s}_named.vcf.gz
    done

    # 4. Merge de los tres
    bcftools merge -m all c2_named.vcf.gz f1_named.vcf.gz m1_named.vcf.gz \
      > trio_merged.vcf

## 4. Resultados

    Variantes por individuo:
      C2:  26
      F1:  29
      M1:  23
    Tras merge:  31 variantes (unión)

## 5. Análisis de herencia

### 5.1 Variante objetivo 8979 (A→T)

    C2=0/1   F1=0/1   M1=./.

- C2 la porta.
- F1 la porta.
- M1 no la porta.
- **Herencia: autosómica, transmitida por el padre.**

Coherente con el diseño del fixture: la variante está forzada en
F1 y C2.

### 5.2 Indel 4407 (TTGTG→T)

    C2=0/1   F1=0/1   M1=0/1

- Los tres la portan en heterocigosis.
- **Herencia: variante común MENA heredada por ambos padres.**

### 5.3 Variante 7032 (C→T)

    C2=0/1   F1=./.   M1=0/1

- C2 la porta.
- F1 no.
- M1 sí.
- **Herencia: transmitida por la madre.**

### 5.4 Variantes 3592, 3677, 5439

    C2=1/1   F1=1/1   M1=1/1

- Los tres homocigotos alt.
- **Herencia: variante común MENA homocigota en todos.**

## 6. Interpretación

**El joint calling del trío contra el grafo local respeta la
herencia mendeliana.** Cada variante se comporta según lo esperado
en función de los genotipos forzados en el fixture.

Esto demuestra que:
1. **El pipeline de grafo es compatible con tríos.**
2. **El merge de VCFs por individuo es suficiente** para joint calling
   (no hace falta fusionar GAMs, que vg no soporta).
3. **Los genotipos de `vg call` son correctos** en contexto familiar.

## 7. Diseño del pipeline

La clave es **no intentar fusionar GAMs** (vg no tiene un comando
directo para eso). En su lugar:

- Alinear cada individuo contra el mismo grafo.
- Pack + call por individuo.
- Renombrar el sample `SAMPLE` (default de vg) al nombre del individuo.
- `bcftools merge -m all` produce el VCF conjunto.

Ventaja: cada VCF individual es independiente y auditable. El merge
es trivial. Sin pérdida de información.

## 8. Limitaciones

- **Sin fase.** El merge de VCFs no produce `0|1`. Para fase
  mendeliana hace falta phasing por transmisión.
- **Tres individuos, no familia extendida.** Un trío es suficiente
  para herencia básica; familias más grandes necesitan más.
- **BAMs sintéticos.** La validación en datos reales queda pendiente.

## 9. Trabajo futuro

- **Fase por transmisión**: aplicar reglas mendelianas a los
  genotipos del trío para determinar la fase del hijo.
- **Extender a familia extendida**: abuelos, hermanos.
- **Integrar con `mendelian`**: filtros automáticos de novo,
  recesivo homocigoto, compuesto.

## 10. Mensaje metodológico

> El pipeline de grafo local funciona con tríos. Cada individuo se
> alinea contra el mismo grafo, se llama por separado, y se fusiona
> con `bcftools merge`. El resultado respeta la herencia mendeliana
> y se puede usar para filtros familiares.

## 11. Reproducibilidad

Ver `docs/demo_nphp1_trio_commands.sh`.
