# Auditoría y formación

Este documento describe cómo se registran las sesiones de
PublicGenomicAgent y cómo se usan para revisión experta y formación
de personal.

Es la base técnica de la idea "chat registrado y auditable" que
sustenta el uso clínico del sistema.

---

## 1. Motivación

Un pipeline genómico produce un VCF. Este sistema produce un caso:
contexto, proceso, decisión, trazabilidad. La trazabilidad no es un
lujo de auditoría — es la diferencia entre una herramienta que da
respuestas y un proceso que se puede revisar.

Tres usos concretos:

1. Revisión de residentes. Un genetista senior revisa 10 casos
   semanales de residentes con transcripts. Ve decisiones, no solo
   diagnósticos.
2. Segunda opinión entre centros. Un centro sin genetista senior
   exporta el transcript. Un centro de referencia lo revisa sin
   volver a procesar los BAMs.
3. Material docente anonimizado. Transcripts anonimizados se
   publican como casos de estudio.

---

## 2. Estructura de una sesión

Cada sesión produce un directorio:

    sessions/<session_id>/
    ├── session.json           # metadatos y eventos
    ├── transcript.md          # versión legible
    ├── inputs/                # referencias a BAMs, VCFs, etc.
    ├── outputs/               # VCFs, reports, HTML
    ├── hashes.txt             # SHA-256 de cada fichero
    └── versions.json          # versiones de tools y binarios

Nunca se copian los BAMs. Solo se referencian por ruta y hash. Los
BAMs originales permanecen donde el usuario los tenga.

---

## 3. Formato de session.json

    {
      "session_id": "2026-09-20T18-14-33_C2-NPHP1",
      "created_at": "2026-09-20T18:14:33Z",
      "user": "genetista_id_x",
      "case_id": "NPHP1-001",
      "consent": {
        "clinical_text_use": true,
        "teaching_use": false,
        "retention_days": 365
      },
      "inputs": [
        {"path": "/data/c2.bam", "sha256": "abc...", "sample": "C2"}
      ],
      "versions": {
        "pga_core": "0.0.1",
        "samtools": "1.19.2",
        "bcftools": "1.19",
        "vg": "1.76.1"
      },
      "events": [
        {
          "ts": "2026-09-20T18:14:35Z",
          "type": "tool_invoked",
          "tool": "fetch_roi",
          "args": {"bam": "c2.bam", "region": "chr2_roi:1-25000"},
          "result": "ok",
          "duration_ms": 342
        }
      ]
    }

---

## 4. Tipos de evento

| Tipo                            | Descripción                                              |
|---------------------------------|----------------------------------------------------------|
| session_start                   | Apertura con case_id y consentimiento                    |
| tool_invoked                    | Ejecución de una tool (args + resultado + duración)      |
| phenotype_extraction_proposed   | Propuesta automática de términos HPO                     |
| phenotype_extraction_confirmed  | Confirmación humana con diff                             |
| mendelian_filter_applied        | Aplicación de regla de segregación                       |
| variant_prioritized             | El genetista marca una variante como prioritaria         |
| variant_filtered_out            | El genetista descarta una variante con motivo            |
| hypothesis_added                | Se añade una hipótesis al caso                           |
| session_end                     | Cierre con resumen                                       |

---

## 5. Transcript legible

transcript.md es una versión humana de session.json. Ejemplo:

    # Sesión 2026-09-20T18-14-33 — NPHP1-001

    Usuario: genetista_id_x
    Caso: NPHP1-001 (trío: C2, F1, M1)
    Consentimiento: uso clínico, sin uso docente, retención 365 días

    ## Contexto

    Trío con fallo renal y enanismo. BAMs alineados a hg38.

    ## Pasos ejecutados

    ### 18:14:35 — fetch_roi (C2, chr2_roi:1-25000) OK 342ms
    ### 18:14:38 — phenotype_extraction_proposed (C2)
      Sistema propone: HP:0000113, HP:0000083
    ### 18:15:02 — phenotype_extraction_confirmed (C2)
      Genetista confirma: HP:0000113
      Corrección: eliminado HP:0000083 (no aplicable)
    ### 18:15:45 — user_annotation
      "descarta variante 7032 por homopolímero"

    ## Resultado final

    Variante candidata: NPHP1 2-110008979-A-T
      GT C2: 0/1  F1: 0/1  M1: ./.
      Herencia: autosómica, transmitida por el padre
      Ranking fenotípico: NPHP1 en posición 1

---

## 6. Anonimización

pga session anonymize <session_id> produce una copia sin:

- Identificadores de paciente (nombres, IDs clínicos).
- Identificadores de genetista.
- Fechas exactas (redondeadas al mes).
- Rutas absolutas (reemplazadas por placeholders).

Mantiene:

- Estructura del caso.
- Decisiones y correcciones.
- Términos HPO.
- Genotipos y variantes.
- Trazabilidad técnica.

El resultado es publicable como caso de estudio sin consentimiento
individual adicional, siempre que el consentimiento original contemple
uso docente anonimizado.

---

## 7. Revisión experta

Un genetista senior puede:

- Anotar una sesión: añadir comentarios en eventos concretos.
- Marcar decisiones como correctas, cuestionables o erróneas.
- Bifurcar el caso: "¿qué pasaría si hubiéramos filtrado X?"
- Exportar como material docente.

La bifurcación es potente: ejecuta la sesión con un filtro diferente
y produce un transcript paralelo. Permite enseñar alternativas.

---

## 8. Reproducibilidad

Dado un session.json y los ficheros de input originales, cualquier
persona puede reproducir la sesión bit a bit:

1. Reconstruir los entornos de micromamba desde los manifests
   (versiones exactas en versions.json).
2. Ejecutar los eventos en orden.
3. Comparar los hashes de los outputs.

Esto convierte cada sesión en un artefacto reproducible. No solo
"esto pasó" — "esto pasa siempre que se ejecute igual".

---

## 9. Retención y derecho al olvido

Cada sesión tiene un retention_days definido en el consentimiento.
Al expirar:

- El session.json y el transcript.md se eliminan.
- Los outputs se eliminan.
- Los hashes se conservan (para auditoría de qué existió).
- Los inputs del usuario no se tocan (son suyos).

Si un paciente retira el consentimiento antes de la expiración:

- Mismo procedimiento, inmediato.
- Se registra el evento de retirada en un log separado.

---

## 10. Cifrado y acceso

En despliegue clínico:

- Los ficheros de sesión se cifran en reposo.
- El acceso se registra (quién lee qué sesión).
- Los transcripts no se exponen en interfaces sin autenticación.

En despliegue de investigación (uso actual):

- Sin cifrado por defecto.
- Acceso local al directorio de sesiones.
- Documentar claramente que no es entorno clínico.

---

## 11. Estado actual

Pendiente de implementación. Este documento fija el formato y la
política, para que cuando se implemente la capa de sesión (tras
mendelian_segregation y phenotype_layer), el diseño esté cerrado.

Componentes previstos:

- pga session start/end
- pga session record (registro automático desde las tools)
- pga session export --format html|pdf|markdown
- pga session anonymize
- pga session diff (comparar dos sesiones)
- pga session replay (reproducir una sesión)

