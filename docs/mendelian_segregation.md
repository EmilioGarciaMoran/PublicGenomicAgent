# Reformulación de la genética mendeliana: de test estadístico a lógica de grafos

Este documento fija el marco conceptual sobre el que se apoya el paso
mendelian del pipeline (local_pangenome.md, sección 8) y la capa de
fenotipo (phenotype_layer.md). Existe porque la intuición de bachillerato
sobre "genética mendeliana" (guisantes de Mendel: muchos individuos, muchas
generaciones, fenotipos discretos, ratios 3:1) no es el modelo correcto para
el problema clínico real: un trío o una familia pequeña, dos o tres
generaciones, fenotipo impreciso, genotipo abierto.

---

## 1. El problema de partida

Un genetista clínico no se sienta con un paciente "para tener un VCF". Se
sienta porque hay un niño con un problema, y su familia está dispuesta a dar
muestras. El razonamiento que sigue no es un experimento controlado con
potencia estadística — es una investigación con pocos datos, dirigida por
hipótesis. Cualquier modelo computacional de este proceso tiene que empezar
por admitir esa diferencia, no camuflarla.

---

## 2. Tres paradigmas — solo uno aplica aquí

### 2.1 Test de asociación poblacional (GWAS)

Compara frecuencia alélica entre casos y controles no emparentados, a
escala de cientos o miles de individuos, y produce un odds ratio con
p-valor. Pregunta: ¿esta variante es más frecuente en enfermos que en
sanos, en la población? No aplica a un trío: no hay "n" con sentido
estadístico.

### 2.2 Análisis de ligamiento clásico (LOD score)

Opera sobre pedigrís, pero necesita muchas meiosis informativas —
históricamente, familias grandes con múltiples afectados a lo largo de
varias generaciones — para tener poder resolutivo. Con un trío (una sola
meiosis por progenitor) el LOD score no aporta nada. Este es el paradigma
que la genética de bachillerato hereda conceptualmente sin decirlo, y es el
paradigma equivocado para el caso clínico habitual.

### 2.3 Segregación mendeliana como consulta lógica

No hay test estadístico sobre la familia. Hay una pregunta binaria, por
variante: ¿es esta configuración de genotipos en padre/madre/hijo
lógicamente consistente con este modelo de herencia? Esto es satisfacción
de restricciones sobre un grafo pequeño, no inferencia estadística. Este
es el paradigma correcto para el problema clínico de trío/familia pequeña.

Precedente directo: GEMINI (Paila et al. 2013, PLoS Comput Biol) trata
un VCF + pedigrí como base de datos relacional y expone auto_dom,
auto_rec, de_novo, comp_het como consultas, no como tests
estadísticos. La herramienta en sí lleva años sin desarrollo activo, pero el
modelo conceptual es exactamente el correcto.

Herramientas de referencia en el ecosistema:

- bcftools +trio-dnm (plugin oficial de bcftools) para detección de
  variantes de novo en tríos.
- GATK PhaseByTransmission para phasing por transmisión.
- WhatsHap --ped para phasing con pedigree.

La implementación de este proyecto es propia y auditable, pero reconoce
estas herramientas como precedentes del dominio.

---

## 3. El pedigrí como grafo dirigido

- Nodos = individuos. Atributos: sexo, estado afectado (booleano
  grueso — ver phenotype_layer.md sección 1 para la distinción con el
  fenotipo fino), genotipo por variante candidata.
- Aristas dirigidas: padre -> hijo, madre -> hijo.
- Librería recomendada: networkx. Un pedigrí de trío o familia pequeña es
  trivial en tamaño; no hay necesidad de una base de datos de grafos
  dedicada.

Ejemplo mínimo:

    import networkx as nx

    pedigree = nx.DiGraph()
    pedigree.add_node("padre", sexo="M", afectado=False)
    pedigree.add_node("madre", sexo="F", afectado=False)
    pedigree.add_node("C2", sexo="M", afectado=True)
    pedigree.add_edge("padre", "C2", rol="padre")
    pedigree.add_edge("madre", "C2", rol="madre")

---

## 4. Reglas de segregación como predicados booleanos

Cada regla es una función (variante, genotipos_por_individuo, pedigree) ->
bool. No hay coeficientes, no hay p-valores.

| Modelo                | Condición                                                              |
|-----------------------|------------------------------------------------------------------------|
| De novo               | Ausente en ambos padres, presente en el hijo afectado                  |
| Recesivo homocigoto   | Heterocigoto en ambos padres, homocigoto alt en el hijo afectado       |
| Recesivo compuesto    | Dos variantes distintas en el mismo gen, cada una de un progenitor     |
| Dominante heredado    | Presente en un progenitor afectado y en el hijo afectado               |
| Ligado a X recesivo   | Hemicigoto en hijo varón afectado, heterocigoto en la madre            |
| Ligado a X dominante  | Presente en ambos sexos, transmisible por ambos progenitores           |
| Mitocondrial          | Transmisión exclusiva por la madre; afectados en ambos sexos           |
| Imprinting            | Dependiente del origen parental. Documentado, no implementado          |

Ejemplo:

    def es_de_novo(genotipo_padre, genotipo_madre, genotipo_hijo):
        return (genotipo_padre == "0/0"
                and genotipo_madre == "0/0"
                and genotipo_hijo in ("0/1", "1/1"))

La función segregation_filter(pedigree, genotipos, modelo) aplica estas
reglas sobre el conjunto de variantes candidatas y devuelve la lista que
las satisface. Es determinista: mismos genotipos, mismo pedigrí, mismo
resultado siempre — coherente con el principio de determinismo de
local_pangenome.md.

---

## 5. Por qué tampoco es "reglas JSON" ingenuas

La lógica de segregación reduce drásticamente el espacio de variantes, pero
no lo suficiente por sí sola: un exoma de trío parte de decenas de miles de
variantes candidatas incluso tras filtros de calidad. El embudo real
combina, en cascada, herramientas de naturaleza distinta:

1. Filtro de calidad (duro) — DP, QUAL, sesgo de hebra.
2. Filtro de frecuencia poblacional (umbral contra referencia externa,
   idealmente la cohorte local de local_pangenome.md en vez de gnomAD
   genérico) — no es un test sobre la familia, es un umbral.
3. Puntuación de consecuencia funcional (VEP + scores tipo CADD/REVEL)
   — modelos estadísticos/ML que actúan como plausibilidad, no como test.
4. Segregación mendeliana (este documento) — lógica dura.
5. Similitud fenotípica (phenotype_layer.md) — ranking blando, no filtro.

Ningún paso aislado resuelve el problema. La combinación en cascada, con
cada paso usando la herramienta adecuada a su naturaleza, es lo que lo
resuelve.

---

## 6. Compound het: el caso que exige más cuidado

Detectar recesivo compuesto correctamente requiere saber que las dos
variantes están en fase distinta (una de cada progenitor), no solo que
ambas existen en el gen. Con datos de trío de lectura corta, esto se puede
inferir en la mayoría de los casos sin necesidad de phasing físico: si una
variante solo aparece en el padre y la otra solo en la madre, y ambas
aparecen en el hijo, la fase queda determinada por herencia biparental. El
caso ambiguo es cuando ambas variantes podrían venir del mismo progenitor
(ej. ambos padres son heterocigotos para ambas posiciones) — ahí la
inferencia por trío no basta y hace falta phasing real (lectura larga, o
métodos estadísticos de read-backed phasing). Esto debe documentarse como
limitación conocida, no resolverse silenciosamente con un supuesto.

Mitigación parcial ya documentada: phasing por trío extendido con
abuelos (local_pangenome.md, sección 13). Con tres generaciones
disponibles, el caso ambiguo de dos padres 0/1 en la misma variante se
resuelve en la mayoría de los casos.

---

## 7. Casos borde a documentar, no a ignorar

- Mosaicismo parental: una variante presente en el hijo a VAF alto
  puede originarse en un progenitor mosaico a VAF bajo, indetectable por un
  variant caller estándar. Esto puede hacer que una variante real parezca
  "de novo" cuando no lo es en sentido estricto. Limitación conocida.

- No-paternidad o no-maternidad biológica: puede aparecer
  incidentalmente al aplicar estas reglas (un patrón de herencia
  imposible bajo el pedigrí declarado). Es una cuestión tanto técnica como
  ética — el sistema debe señalar la inconsistencia sin inferir ni mostrar
  conclusiones sobre paternidad/maternidad de forma automática o expuesta
  sin mediación clínica.

---

## 8. Integración con el resto del pipeline

    fetch_roi -> qc_bam -> local_pangenome -> call_variants
                                                  |
                                                  v
                                       mendelian_segregation (este documento)
                                                  |
                                       lista corta de genes candidatos
                                                  |
                                                  v
                                          phenotype_layer.md
                                                  |
                                  evidencia combinada, presentada, no decidida

---

## 9. Validación

Sobre el fixture NPHP1 (local_pangenome.md, sección 9):

- test_segregation_detects_de_novo
- test_segregation_detects_auto_rec_homocigoto
- test_segregation_detects_comp_het_fase_biparental
- test_segregation_flags_comp_het_fase_ambigua_como_no_resuelto
- test_segregation_flags_pedigree_inconsistente
- test_segregation_requires_prior_quality_and_freq_filters

---

## 10. Decisiones de diseño abiertas

- Cómo representar el estado "no resuelto por trío, requiere phasing
  físico" en la salida de compound_het.
- Umbral de VAF para sospechar mosaicismo parental no detectado.
- Protocolo explícito de manejo ante inconsistencia de pedigrí detectada
  (a quién se notifica, cómo se documenta, qué no se muestra
  automáticamente): requiere revisión por genetista clínico.
