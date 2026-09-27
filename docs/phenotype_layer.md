# Capa de fenotipo

Este documento describe el diseño de la capa de fenotipo de
PublicGenomicAgent. Es complementaria a local_pangenome.md y a la
lógica de segregación mendeliana (mendelian_segregation.md): mientras
esas dos capas trabajan sobre genotipo y pedigrí, esta capa incorpora
la evidencia clínica del paciente para priorizar, nunca para filtrar
de forma dura, la lista de genes candidatos.

---

## 1. Motivación

Un caso clínico real (trío, fenotipo impreciso, genotipo abierto) no
tiene potencia estadística suficiente para tratar el fenotipo como una
variable poblacional. Pero descartar el fenotipo por completo
desperdicia la señal más rica que aporta el clínico: la descripción
del problema del paciente.

El error habitual es tratar el fenotipo como una sola cosa. En
realidad hay dos fenotipos distintos, con dos usos distintos:

- Fenotipo grueso (afectado / no afectado): un atributo binario por
  individuo, que ya vive en el grafo de pedigrí y alimenta la lógica
  de segregación (mendelian_segregation.md). No requiere esta capa.
- Fenotipo fino (términos HPO — "insuficiencia renal", "poliuria",
  "quistes corticomedulares"): una descripción semántica rica, por
  individuo, que no filtra variantes, sino que prioriza genes ya
  reducidos por la segregación, comparándolos contra una base de
  conocimiento externa (HPO -> gen -> enfermedad).

Esta capa se ocupa exclusivamente del segundo caso.

---

## 2. Filosofía

- Extracción auditable, no conversacional: la interpretación de un
  relato clínico es un artefacto versionado y confirmable, nunca un
  paso oculto dentro de una respuesta de chat.

- El fenotipo prioriza, no filtra: nunca se usa como exclusión dura
  de una variante, salvo que el propio clínico decida excluir un gen
  explícitamente (acción registrada como decisión humana, no del
  sistema).

- La negación es una señal, no ruido: "sin anomalías cardíacas" no es
  lo mismo que ausencia de mención.

- Consistencia entre afectados como evidencia de primera clase: si
  hay más de un afectado en la familia, comparar sus perfiles
  fenotípicos contra el mismo gen candidato es una señal gratuita que
  la mayoría de pipelines no explota.

- Sin dependencia de nube: extractores y motores de ranking corren en
  local.

- El sistema no decide, presenta evidencia.

- Consentimiento explícito: el uso del relato clínico para extracción
  HPO requiere consentimiento informado explícito, registrado en el
  session_log. Sin consentimiento, la capa no se activa.

---

## 3. Entradas

- Relato clínico por individuo afectado: texto libre, potencialmente
  multilingüe (árabe, español, inglés).
- Alternativamente, términos HPO ya codificados si el clínico los
  introduce directamente.
- Pedigrí con estado afectado/no afectado.
- Lista corta de genes candidatos producida por la capa de
  segregación mendeliana (entrada obligatoria; esta capa no opera
  sobre el genoma completo).

---

## 4. Salidas

- hpo_terms/<individuo>.json: términos HPO confirmados por
  individuo, con proveniencia completa:

    {
      "individuo_id": "C2",
      "terms": [
        {"hpo_id": "HP:0000113", "label": "Enfermedad renal poliquística",
         "negated": false, "confidence": 0.87, "source_span": [42, 68],
         "extractor": "llm:qwen3.5:9b", "confirmed_by": "genetista_id_x"}
      ]
    }

- gene_ranking/<individuo>.tsv: genes candidatos rankeados por
  similitud fenotípica para ese individuo (salida de Phen2Gene).

  Formato: gene, rank, score, matched_hpo_terms (tab-separated).

- consistency_report.json: comparación del ranking entre todos los
  afectados de la familia — ver sección 7.

- Eventos en session_log (ver audit_and_training.md):
  phenotype_extraction_proposed, phenotype_extraction_confirmed.

---

## 5. Herramientas (verificadas)

| Función                     | Herramienta              | Offline | Instalación      |
|-----------------------------|--------------------------|---------|------------------|
| Texto -> HPO (inglés)       | txt2hpo                  | Sí      | pip, CPU         |
| Texto -> HPO (inglés, mejor)| pleio-hpo                | Sí      | pip, ~1.5 GB     |
| Texto -> HPO (multilingüe)  | LLM local + confirmación | Sí      | ya integrado     |
| HPO -> genes candidatos     | Phen2Gene                | Sí      | conda            |

Evitado explícitamente:

- Doc2HPO estándar usa una API web pública — el texto clínico sale
  de la máquina. Solo aceptable con instalación local propia y
  licencia UMLS.
- PhenoGPT: peor F1 del benchmark (0.341) con la mayor huella
  (~22 GB, GPU). No apto para portátil.

No existe extractor dedicado maduro para árabe o español. Para esos
idiomas, LLM local con confirmación humana obligatoria.

---

## 6. Algoritmo

Paso 1. Ingesta del relato clínico, uno por individuo afectado.

Paso 2. Extracción propuesta:

- Texto en inglés: txt2hpo o pleio-hpo (deterministas, versión
  fija, sin LLM).
- Otro idioma: LLM local, con prompt que exige detección explícita
  de negación y devuelve el span de texto origen de cada término.

Paso 3. Confirmación humana obligatoria. El genetista ve la
propuesta, añade, elimina o corrige términos. Nada pasa a los pasos
siguientes sin confirmación. Se registra el diff completo
(terms_proposed vs terms_final) en session_log.

Paso 4. Ejecutar Phen2Gene con los términos confirmados de cada
afectado por separado. No se fusionan los perfiles.

Paso 5. Cruzar cada gene_ranking/<individuo>.tsv con la lista corta
de genes candidatos por segregación (intersección, no unión).

Paso 6. Generar consistency_report.json (sección 7).

Paso 7. Registrar todos los artefactos y decisiones en session_log.

---

## 7. Consistencia entre afectados

Con más de un afectado en la familia, cada uno con fenotipo
ligeramente distinto, se puede comprobar si un gen candidato explica
consistentemente a todos los afectados, no solo al probando.

    {
      "gene": "NPHP1",
      "ranking_por_individuo": {
        "C2": {"posicion": 1, "score": 0.91},
        "C3": {"posicion": 14, "score": 0.22}
      },
      "alerta": "divergencia_alta",
      "interpretacion_sugerida": "posible heterogeneidad de locus o hallazgo incidental en C2 — revisar manualmente"
    }

Esta alerta no es una conclusión automática — es evidencia que se
presenta al genetista.

---

## 8. Manejo de negación

Un relato clínico del tipo "sin anomalías cardíacas, sin hipoacusia"
invierte completamente el significado si un extractor ingenuo marca
esos términos como presentes. txt2hpo y pleio-hpo manejan negación
sobre inglés (arquitecturas tipo NegEx). Para extracción vía LLM, el
prompt debe exigir el campo negated, y el paso de confirmación humana
es la última red de seguridad.

---

## 9. Limitación de idioma

No hay extractor dedicado maduro para árabe o español. Para esos
idiomas, la capa depende del LLM (no determinista). Consecuencia: la
confirmación humana no es opcional.

Además, dentro del árabe hay variación dialectal (golf, levantino,
magrebí) que puede afectar la extracción. La detección de idioma y
dialecto es un paso previo recomendable.

---

## 10. Esquema de datos

    individuals
    ├─ id
    ├─ pedigree_id
    ├─ afectado (bool)
    └─ hpo_terms[]

    gene_phenotype_kb (interno a Phen2Gene, no se reimplementa)
    ├─ gen
    ├─ hpo_term_id
    └─ enfermedad (OMIM/Orphanet)

No se crea una tabla gene_phenotype_kb propia: Phen2Gene ya encapsula
phenotype_to_genes.txt internamente.

---

## 11. Integración con el resto del pipeline

    mendelian_segregation.md  -> lista corta de genes candidatos (lógica dura)
    phenotype_layer.md        -> ranking de esos mismos genes (similitud fenotípica)
                                  + consistency_report entre afectados
                               -> evidencia combinada, presentada, no decidida

Materialización concreta: el ranking fenotípico se presenta como un
orden de genes candidatos dentro de la lista corta producida por
segregación. Cada gen lleva asociado:

- Posición en el ranking para cada afectado.
- Score de similitud para cada afectado.
- Alerta de divergencia si los afectados discrepan.
- Términos HPO que soportan el ranking.

No se convierte en un score único por variante.

---

## 12. Validación

Extensión del fixture NPHP1:

- test_phenotype_extraction_handles_negation
- test_phen2gene_ranks_nphp1_high
- test_consistency_report_flags_divergence
- test_no_extraction_without_confirmation
- test_no_extraction_without_consent

---

## 13. Decisiones de diseño

- Versión de HPO fijada: v2025-XX-XX (fecha exacta en envs_pin.lock).
  Actualización trimestral con test de regresión.
- Extractor por defecto según idioma: pendiente de definir umbral de
  confianza en la detección de idioma.
- Representación de "no evaluado" frente a "ausente confirmado": HPO
  distingue excluded de ausencia de mención. Se modela como tercer
  estado explícito (not_assessed).
- Fusión de cohortes de fenotipo (Al Mena, EGP, gnomAD MID):
  pendiente, en paralelo a local_pangenome.md sección 10.

---

## 14. Requisitos de hardware (actualizado)

La extracción de fenotipos con RDMA **no cabe en un MacBook Air de 16 GB**.

### Números reales

El modelo recomendado por RDMA es **Mistral Small 24B en 4-bit**:

| Componente | VRAM/RAM |
|---|---|
| Pesos del modelo (24B, 4-bit NF4) | ~14 GB |
| KV cache (contexto 8k) | ~1 GB |
| Overhead del runtime (PyTorch, bitsandbytes) | ~2-3 GB |
| **Total mínimo** | **~17-18 GB** |

Esto descarta:
- MacBook Air 16 GB RAM unificada.
- Cualquier portátil con GPU < 20 GB VRAM.

### Hardware mínimo viable

- **GPU NVIDIA RTX 3090 (24 GB VRAM)** — recomendada explícitamente por el paper de RDMA.
- **GPU NVIDIA RTX 4090 / A5000 (24 GB VRAM)** — equivalente.
- **Mac con ≥32 GB RAM unificada** — vía backend `llama_cpp` con GGUF (no probado aún).
- **Servidor con A100/H100** — para uso hospitalario.

### Modelos más pequeños: no es una opción

El paper de RDMA es explícito: **por debajo de 24B, los modelos no siguen el "structured output" de forma fiable** y la precisión cae drásticamente. No es aceptable para uso clínico.

### Consecuencia arquitectónica

El sistema PublicGenomicAgent tiene dos naturalezas:

**Capa de variantes (BAMs, ROIs, grafos)**:
- Cabe en un portátil modesto (16 GB RAM).
- Funciona sin GPU.
- `fetch_roi`, `qc_bam`, `call_variants`, `local_pangenome`, `build_local_graph`, `align_to_graph`, `call_from_graph`, `compare_vcfs`, `mendelian_filter`.

**Capa de fenotipo (extracción HPO/ORPHA)**:
- Requiere GPU ≥20 GB VRAM **o** servidor dedicado.
- No cabe en un portátil de 16 GB.
- `extract_hpo` (RDMA).

Para un despliegue hospitalario (KFSHRC, KAUST), la arquitectura natural es:

    [MacBook Air del genetista]                    [Servidor GPU]
         │                                              │
         │  fetch_roi, vg, call_variants                │
         │  (todo lo que no necesita LLM)               │
         │                                              │
         │  texto clínico ──────────────────────────────>
         │                                              │  RDMA
         │                                              │  (Mistral 24B 4-bit)
         │  <──────────────────── hpo_terms.json ───────│
         │                                              │
         ▼                                              ▼

El Mac es el terminal clínico. El servidor GPU es el motor de extracción
de fenotipos. Los datos nunca salen de la red local.

Esta división es coherente con los requisitos de privacidad del contexto
saudí: los datos clínicos de pacientes no pueden salir del país ni
enviarse a APIs externas.

### Estado actual del código

`extract_hpo` (tool `pga tool extract-hpo`) soporta 5 backends:

- `local` → HuggingFace + bitsandbytes. Requiere GPU ≥20 GB.
- `llama_cpp` → GGUF local. Requiere ~14 GB RAM. No probado aún.
- `openrouter` → API externa. **No apto para datos clínicos reales.**
- `api` (Groq) → API externa. **No apto para datos clínicos reales.**
- `azure` → API externa. **No apto para datos clínicos reales.**

Los tres últimos son válidos **solo para pruebas y desarrollo**.

### Verificación sin GPU

Sin GPU, `extract_hpo` falla en **menos de 10 segundos** con un mensaje
claro:

    RuntimeError: Backend 'local' requiere GPU NVIDIA con CUDA.
    No se detectó CUDA disponible. Alternativas: --backend openrouter
    o --backend api, o usar una máquina con GPU (>=20 GB VRAM para
    Mistral 24B).

Esto es deliberado: el sistema **no intenta** cargar 14 GB de modelo en
CPU ni descargar 48 GB de pesos. Falla antes.

### Trabajo futuro (arquitectura de servidor)

Para el despliegue hospitalario real, las piezas que faltarán son:

1. **Servidor RDMA**: servicio FastAPI que expone `POST /extract_hpo`
   con el texto clínico y devuelve el JSON de términos HPO.
2. **Backend `remote` en `extract_hpo`**: cliente HTTP que envía el
   texto al servidor y recibe la respuesta.
3. **Autenticación y cifrado**: TLS + tokens, para red hospitalaria.
4. **Docker/compose**: empaquetado del servidor para despliegue en
   KFSHRC, KAUST, etc.

Estas piezas se implementarán cuando haya un caso de uso real con
hardware disponible.
