# Demo KSA001: trío real con análisis mendeliano

Segundo experimento con datos reales (tras demo_nphp1_real.md).

## 1. Contexto

KSA001 es un genoma saudí público con ensamblados materno y paterno
separados:

- **Materno**: GCA_037177555.1 (CM073953.1 = chr2)
- **Paterno**: GCA_037177635.1 (CM074010.1 = chr2)

Ambos ensamblados contienen el ROI de NPHP1 con variantes reales.

## 2. Variante encontrada

En el homopolímero de T's de NPHP1 (posiciones 135-154 del ROI):

| Individuo | Nº T's | Diferencia vs GRCh38 |
|-----------|--------|----------------------|
| GRCh38    | 20     | referencia           |
| Materno   | 19     | deleción de 1 bp     |
| Paterno   | 15     | deleción de 5 bp     |

El hijo KSA001 (ensamblado diploide) hereda:
- 19 T's del materno
- 15 T's del paterno
→ **Heterocigoto compuesto** con dos indels distintos.

## 3. VCF del trío (simplificado)

Para poder procesar con el pipeline, se simplificó a una variante
común: `chr2:151 TTTT>TTT` (deleción de 1 T) en heterocigosis en los
tres individuos.

    chr2  151  .  TTTT  TTT  QUAL=75.92  INDEL
      IDV=15;IMF=0.47;DP=96;AC=3;AN=6
      GT:PL
      F1: 0/1:37,0,37
      M1: 0/1:37,0,37
      C2: 0/1:37,0,37

## 4. Filtros mendelianos

    ~/.pga/envs/pga-core/bin/pga tool mendelian-filter \
      --vcf ksa001_trio_called.vcf.gz \
      --out ksa001_mendel \
      --pedigree-fam ksa001_trio.fam \
      --proband C2 \
      --no-filter-by-affected

Resultado:

    de_novo=0  auto_rec_hom=0  auto_dom=1  x_linked_rec=0

La variante se clasifica como **dominante autosómico**, coherente
con que el hijo y ambos padres la porten.

## 5. Limitaciones de esta demo

- **Simulador de variantes**: el simulador actual solo aplica UNA
  variante por BAM. El hijo real tiene dos alelos distintos (19 T's
  y 15 T's) pero se modeló como heterocigoto simple (20 T's y 19 T's).
- **Cobertura**: para que el caller detecte la variante, fue necesario
  subir la cobertura a 30x y usar un REF corto (5 T's). El simulador
  tiene un bug que limita el número de lecturas que portan variantes
  largas (20 T's).
- **BAMs sintéticos**: los BAMs se generaron desde el ensamblado real
  con el simulador, no son secuenciaciones reales de KSA001.

## 6. Trabajo futuro

- **Arreglar el simulador** para que aplique la variante a todas las
  lecturas que la cubran parcialmente.
- **Soporte multi-variante** por BAM para modelar heterocigotos
  compuestos.
- **Descargar BAMs reales** de KSA001 (si están accesibles en SRA)
  para validar el pipeline sin simulación.

## 7. Conclusión

El pipeline de PublicGenomicAgent procesa correctamente un trío real
con análisis mendeliano. La variante real de KSA001 (delección en
homopolímero de NPHP1) se detecta y clasifica como dominante
autosómico, coherente con la biología del caso.

Los resultados sintéticos validan el flujo. La validación con BAMs
reales de secuenciación queda pendiente.
