# PublicGenomicAgent — pitch

**A privacy-first, laptop-scale copilot for clinical genomics.**

Dado un caso clínico (fenotipo, familia, BAMs ya alineados),
PublicGenomicAgent completa el análisis genómico end-to-end en un
portátil de consumo. Sin clúster, sin GPU, sin red. Y cuando el
pipeline lineal falla —duplicaciones segmentarias— activa pangenoma
local dirigido para recuperar las variantes patogénicas que se
pierden.

---

## El problema

1. **El cuello de botella no es el alineamiento, es la interpretación.**
   Los hospitales ya producen BAMs con DRAGEN, BWA+GATK o minimap2.
   Lo que falta no es otra pipeline de alineamiento, sino un copiloto
   que ayude al genetista a pasar de "tengo un BAM" a "tengo una
   decisión clínica" —sin necesidad de infraestructura pesada.

2. **Las duplicaciones segmentarias engañan al pipeline lineal.**
   En regiones con LCRs (NPHP1, CYP2D6, LPA, KIR), el alineamiento
   lineal produce MAPQ=0 y los callers rechazan variantes patogénicas
   reales. La referencia lineal (GRCh38) está sesgada hacia
   poblaciones europeas, penalizando cohortes del Golfo y MENA.

3. **Las herramientas de análisis genómico asumen clúster, GPU y
   conexión a la nube.**
   El genetista clínico no puede esperar 6 horas para un caso, ni
   subir datos de pacientes a un servidor externo. Necesita un
   copiloto local, rápido y auditable.

---

## Los tres pilares

### 1. Análisis genético completo en un portátil

**No alineamos.** Partimos de BAMs ya alineados por el pipeline
habitual del hospital. El valor del sistema está en la capa de
interpretación dirigida por fenotipo, no en el alineamiento.

Dado un caso clínico —fenotipo, familia, BAMs—, PublicGenomicAgent
completa el análisis en un ordenador de consumo:

```
Texto clínico (opcional)
    ↓
Extracción HPO (extract_hpo)
    ↓
Priorización fenotípica (phenotype_ranking, LIRICAL)
    ↓
ROI dirigido (resolve_input → fetch_roi)
    ↓
QC + joint calling (qc_bam, call_variants, joint_call)
    ↓
Filtros mendelianos (mendelian_filter)
    ↓
Anotación ClinVar (annotate_variants)
    ↓
Priorización clínica (prioritize_variants)
    ↓
Informe HTML + IGV
```

**Números**:
- Un caso de trío real (Osteopetrosis, 1 Mb ROI) se procesa en
  **19.6 segundos**, incluyendo:
  - 3 QC de BAMs de ~350 MB
  - 3 sub-BAMs por ROI
  - 3 variant callings individuales
  - 1 joint calling multi-sample
  - 1 filtro mendeliano (de novo, recesivo, dominante, X-linked)
  - 1 anotación ClinVar
  - 1 priorización clínica
- **13 tool calls**, todas trazables, todas deterministas.
- **Sin GPU, sin clúster, sin red.**

**Auditoría**:
- Cada tool call queda registrada en `session.json` con input
  validado, output serializado, duración y error.
- Reproducible: mismo input → mismo output.

---

### 2. Pangenoma local dirigido: recuperar lo que el lineal pierde

**No pangenomizamos todo el genoma.** El sistema decide cuándo el
pangenoma aporta valor. Es una herramienta **focalizada y
condicionada**, no una bala de plata.

**Los dos benchmarks**:

**Benchmark A — trío real sin duplicaciones segmentarias
(Osteopetrosis, chr8)**:

| Pipeline | Variantes detectadas |
|---|---|
| Lineal (`bcftools call`) | 102 |
| Pangenoma local (`vg`) | 91 |
| Compartidas | 76 |
| Concordancia | 83.5% |

**Conclusión**: en regiones sin SDs, los dos pipelines son
comparables. El pangenoma es más conservador (excluye ~20
artefactos de baja calidad que el lineal sí incluye), pero no
añade variantes nuevas.

**Benchmark B — duplicación segmentaria sintética (chr1)**:

Referencia con una duplicación en tándem de 100 bp. BAM del
proband con 50 reads que cubren la región duplicada, todas
llevando la variante `chr1:1050 T>G` (hom alt). Todas las reads
con MAPQ=0 (simulan la ambigüedad de alineamiento).

| Pipeline | Variantes detectadas | Genotipo |
|---|---|---|
| Lineal (default) | 0 | — |
| Lineal (sin filtro MAPQ) | 0 | — |
| **Pangenoma local (vg)** | **1** | **1/1, DP=49, GQ=132** |

**Conclusión**: en regiones con SDs, el pipeline lineal **pierde
la variante patogénica**. El pangenoma local la recupera con
confianza alta. Es una diferencia medida, no teórica.

**El argumento completo**: el pangenoma no es siempre mejor. Es
**decisivo** en las regiones donde importa (SDs, LCRs) y **neutro**
en el resto. El sistema lo sabe y lo activa cuando aporta.

---

### 3. Interfaz conversacional y trazabilidad completa

**La UI es la puerta de entrada. El valor está en la orquestación.**

El usuario describe el caso en lenguaje natural:
> "Tengo un probando de 8 años, hijo de primos hermanos, con
> poliuria, polidipsia e insuficiencia renal progresiva. Aquí
> están los BAMs del trío."

El agente:

1. Extrae términos HPO del texto (`extract_hpo`).
2. Prioriza enfermedades candidatas (`phenotype_ranking`, LIRICAL).
3. Construye el `CaseManifest` (pedigrí, ROIs, proband).
4. Orquesta la secuencia de tools adecuada para el caso.
5. Genera un informe HTML con pedigree, variantes candidatas,
   ranking clínico, anotación ClinVar y enlaces a IGV.
6. Emite un `session.json` con la traza completa.

**Los LLM son locales (Ollama).** Sin red, sin proveedores externos.
Reglas deterministas primero; el LLM solo interviene cuando las
reglas no saben qué hacer. Validación Pydantic de cada decisión.
Fallback determinista si el LLM falla.

**Privacidad por arquitectura**:
- Prompt al LLM solo con basenames (`proband.bam`, no rutas absolutas).
- Ollama en `localhost`, sin llamadas externas.
- `NullLLMClient` para forzar funcionamiento determinista sin LLM.
- Cada garantía está cubierta por un test unitario.

---

## Fábrica de casos: el proyecto Sandbox

PublicGenomicAgent no depende de datos reales para validarse.
Un **proyecto paralelo** (Sandbox) genera catálogos de casos
sintéticos identificados por **SPDI canónico**, agrupados en
tríos, con BAMs reales y ground truth.

**Separación limpia por contrato**:

- Sandbox **genera** los datos (responsabilidad del paralelo).
- PGA **consume** los datos (responsabilidad del core).
- El único punto de acoplamiento es `manifest.yaml`.
- El generador es opaco para PGA: no sabe ni necesita saber
  cómo se produce cada caso.

**Casos actuales en el Sandbox**:

| Gen | SPDI | Fenotipo | Proband |
|---|---|---|---|
| NPHP1 | `NC_000002.12:110800000:290000:` | Nefronoptisis | HOM_DEL |
| CYP2D6 | `NC_000022.11:42128940:C:T` | Farmacogenómica | HOM_ALT |
| CONTROL_NPHP1 | `NC_000002.12:110800000::` | Control sano | HOM_WT |

**Consumo desde PGA**:

```bash
pga sandbox list ~/sandbox
pga sandbox describe ~/sandbox/genes/NPHP1/NC_000002.12_110800000_290000_
pga sandbox run ~/sandbox/genes/CYP2D6/NC_000022.11_42128940_C_T \
    --reference ~/sandbox/test_ref.fa
```

**Por qué importa para el pitch**:

1. **Validación sin datos de pacientes.** No necesitamos acceso
   a cohortes reales para verificar que el pipeline funciona.
2. **Reproducibilidad total.** Cada caso está identificado por
   SPDI canónico. Dos máquinas generan el mismo caso.
3. **Extensibilidad.** Añadir un caso nuevo es crear un
   directorio, no modificar PGA.
4. **Ground truth explícito.** Cada manifest define el
   resultado esperado (`classification`, `genotypes`). Los
   tests comparan contra eso automáticamente.

**Verificado end-to-end**:

```
Caso CYP2D6 → PGA run-trio → 1 variante recesiva detectada
  padre 0/1, madre 0/1, proband 1/1  →  auto_rec_hom  →  score 3.0
```

Cubierto por `tests/integration/test_sandbox_run.py`.

## Fábrica de casos: el proyecto Sandbox

PublicGenomicAgent no depende de datos reales para validarse.
Un **proyecto paralelo** (Sandbox) genera catálogos de casos
sintéticos identificados por **SPDI canónico**, agrupados en
tríos, con BAMs reales y ground truth.

**Separación limpia por contrato**:

- Sandbox **genera** los datos (responsabilidad del paralelo).
- PGA **consume** los datos (responsabilidad del core).
- El único punto de acoplamiento es `manifest.yaml`.
- El generador es opaco para PGA: no sabe ni necesita saber
  cómo se produce cada caso.

**Casos actuales en el Sandbox**:

| Gen | SPDI | Fenotipo | Proband |
|---|---|---|---|
| NPHP1 | `NC_000002.12:110800000:290000:` | Nefronoptisis | HOM_DEL |
| CYP2D6 | `NC_000022.11:42128940:C:T` | Farmacogenómica | HOM_ALT |
| CONTROL_NPHP1 | `NC_000002.12:110800000::` | Control sano | HOM_WT |

**Consumo desde PGA**:

```bash
pga sandbox list ~/sandbox
pga sandbox describe ~/sandbox/genes/NPHP1/NC_000002.12_110800000_290000_
pga sandbox run ~/sandbox/genes/CYP2D6/NC_000022.11_42128940_C_T \
    --reference ~/sandbox/test_ref.fa
```

**Por qué importa para el pitch**:

1. **Validación sin datos de pacientes.** No necesitamos acceso
   a cohortes reales para verificar que el pipeline funciona.
2. **Reproducibilidad total.** Cada caso está identificado por
   SPDI canónico. Dos máquinas generan el mismo caso.
3. **Extensibilidad.** Añadir un caso nuevo es crear un
   directorio, no modificar PGA.
4. **Ground truth explícito.** Cada manifest define el
   resultado esperado (`classification`, `genotypes`). Los
   tests comparan contra eso automáticamente.

**Verificado end-to-end**:

```
Caso CYP2D6 → PGA run-trio → 1 variante recesiva detectada
  padre 0/1, madre 0/1, proband 1/1  →  auto_rec_hom  →  score 3.0
```

Cubierto por `tests/integration/test_sandbox_run.py`.

## Qué nos diferencia

| Aspecto | Pipeline clásico | Plataforma cloud | **PublicGenomicAgent** |
|---|---|---|---|
| Punto de partida | FASTQ / exoma | FASTQ / exome | **BAM + fenotipo** |
| Scope | Todo el genoma | Todo el genoma | **ROI dirigido por clínica** |
| Infraestructura | Clúster + GPU | Cloud + APIs | **Portátil de consumo** |
| Tiempo por caso | Horas | Minutos + upload | **19.6 s (trío real)** |
| Privacidad | Local si hay clúster | Datos fuera | **Air-gapped por diseño** |
| Salida | VCF exhaustivo | VCF + panels | **Decisión clínica + traza** |
| Pangenoma | Opcional | A veces | **Dirigido, cuando aporta** |
| Precio | Coste de infraestructura | Suscripción | **Open source (Apache 2.0)** |

---

## A quién va dirigido

**Genetistas clínicos** que necesitan analizar casos mendelianos
de forma rápida y auditable, sin depender de infraestructura pesada.

**Grupos de investigación** en consanguinidad y enfermedades
recesivas de poblaciones del Golfo y MENA, donde el sesgo de
referencia lineal es un problema real y medible.

**Hospitales con soberanía de datos** (KFSH&RC, KAUST, QGP,
hospitales nacionales) que requieren que ningún dato salga del
entorno controlado.

**Investigadores en pangenómica clínica** que quieren medir el
valor real de los pangenomas locales con benchmarks reproducibles.

---

## Estado del proyecto

**Fase 1 — Consolidación** ✅ completada
- Repositorio público, Apache 2.0, CI en GitHub Actions.
- Documentación: README, quickstart, privacy, known issues, backlog.
- Instalador unificado (`setup.sh`).

**Fase 2 — Valor clínico** ✅ completada
- 16 tools registradas, todas con Pydantic I/O.
- Pipeline completo de 13 pasos (con `--clinvar`).
- Informe HTML con pedigree, variantes, ranking, IGV.
- Integración con ClinVar real (100 MB).
- Priorización clínica con scoring.
- Integración con IGV Desktop (session.xml + control URLs).

**Fase 3 — Datos reales y benchmarks** 🚧 en curso
- ✅ Benchmark lineal vs pangenoma (Osteopetrosis).
- ✅ Benchmark duplicación segmentaria sintética.
- ✅ Primer run sobre BAMs reales (19.6 s).
- ❌ Trío real con cobertura del gen causal (bloqueado por
  descarga de datos).
- ❌ Test de asociación genotipo-fenotipo formal.
- ❌ Pitch formalizado para instituciones saudíes.

**Fase 4 — Investigación** ❌ pendiente
- HPO sin GPU.
- Publicación (Bioinformatics, GigaScience).
- Validación clínica en colaboración con hospitales.

**Tests**: 120 passed (85 unit, 35 integration).

---

## Roadmap para colaboración

**Si eres un hospital / institución (KFSH&RC, KAUST, QGP)**:
- Necesitamos acceso a tríos reales con fenotipo claro.
- Ofrecemos: análisis local, sin exfiltración de datos.
- Colaboración: validación clínica + publicación conjunta.

**Si eres un grupo de investigación**:
- Podemos compartir benchmarks sintéticos reproducibles.
- Podemos colaborar en extender el benchmark a NPHP1, CYP2D6, LPA.
- Podemos publicar la metodología como artifact reproducible.

**Si quieres contribuir código**:
- Repo público en GitHub: `EmilioGarciaMoran/PublicGenomicAgent`.
- Apache 2.0.
- Backlog abierto y priorizado (`docs/backlog.md`).
- CI en cada push.

---

## Lo que NO somos

Es igual de importante aclarar los límites:

- **No alineamos.** Partimos de BAMs ya alineados.
- **No reemplazamos al pipeline clínico.** Lo complementamos.
- **No diagnosticamos.** Priorizamos y trazamos; el diagnóstico
  es del genetista.
- **No usamos LLMs remotos.** Todo local.
- **No somos una bala de plata pangenómica.** El pangenoma se
  activa solo cuando aporta.
- **No soportamos cohortes GWAS.** El enfoque es caso-por-caso,
  no asociación a nivel de población.

---

## Contacto

**Repositorio**: https://github.com/EmilioGarciaMoran/PublicGenomicAgent

**Documentación técnica**:
- [`README.md`](../README.md)
- [`docs/demo_osteo_real.md`](demo_osteo_real.md) — primer caso real.
- [`docs/benchmark_linear_vs_pangenome.md`](benchmark_linear_vs_pangenome.md) — benchmark A.
- [`docs/benchmark_duplications.md`](benchmark_duplications.md) — benchmark B.
- [`docs/privacy.md`](privacy.md) — política de privacidad.
- [`docs/backlog.md`](backlog.md) — estado del proyecto.

---

## Elevator pitch (30 segundos)

> "PublicGenomicAgent es un copiloto genómico clínico para
> portátiles. Dado un caso clínico y BAMs ya alineados, completa
> el análisis end-to-end —sin alinear, sin clúster, sin GPU—:
> desde la descripción del fenotipo hasta el ranking
> genotipo-fenotipo, con anotación ClinVar y trazabilidad completa.
> Cuando el pipeline lineal falla por duplicaciones segmentarias,
> activa pangenoma local dirigido y recupera las variantes
> patogénicas que se pierden. Todo local, sin red, con privacidad
> por arquitectura."

