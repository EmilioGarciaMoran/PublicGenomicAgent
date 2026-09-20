# Casos de uso

Recordatorio de los cuatro escenarios que el sistema debe soportar.
Borrador; desarrollo completo en sesión dedicada.

Todos parten del mismo `CaseManifest`. Se diferencian en qué campos
tienen rellenos, qué tools se invocan y qué falta por construir.

---

## Caso 1 — Genetista clínico con historias/informes

**Input:** PDFs, historias clínicas, notas. A veces con fenotipo
explícito, a veces sin él.

**Necesidad:**
- Interrogar hasta asegurar que hay fenotipo claro.
- Extraer pedigree del texto (preguntar si hace falta).
- Proponer ROIs para encargar un panel de genes.
- Ofrecer conocimiento OMIM/Orphanet como contexto clínico.

**Output esperado:** lista de genes/regiones + contexto clínico +
propuesta de panel de secuenciación.

**Campos del CaseManifest:** pedigree (extraído), hpo_terms
(interrogado), candidate_rois (derivados por HPO → genes).

**Falta por construir:** carga de PDFs, modo interrogación (capa
agéntica), HPO → genes (Phen2Gene), cliente OMIM/Orphanet.

---

## Caso 2 — El mismo genetista vuelve con los BAMs

**Input:** los mismos datos clínicos + BAMs alineados.

**Necesidad:**
- Control de calidad.
- Detección de variantes.
- Control visual (IGV).
- Alternativas contra reference bias (grafo local).

**Output esperado:** VCF anotado + comparación lineal vs grafo +
reporte IGV interactivo.

**Campos del CaseManifest:** pedigree, hpo_terms, candidate_rois.

**Falta por construir:** poco. Las tools ya existen. Falta
integración en un flujo único.

---

## Caso 3 — Laboratorio de secuenciación poblacional

**Input:** BAMs de una cohorte, sin historia clínica.

**Necesidad:**
- Análisis de regiones complejas.
- Verificación de variantes.
- Descubrimiento sin fenotipo.

**Output esperado:** VCFs por muestra + comparaciones, sin capa de
priorización fenotípica.

**Campos del CaseManifest:** candidate_rois (dados), sin pedigree,
sin HPO.

**Falta por construir:** procesamiento en lote (loop sobre BAMs),
merge de VCFs.

---

## Caso 4 — GWAS sobre coordenadas concretas

**Input:** muchos BAMs, coordenadas definidas (no fenotipo).

**Necesidad:** estadística poblacional sobre regiones, comparación
caso/control si aplica.

**Output esperado:** VCF conjunto + análisis PLINK + reporte de
asociación.

**Campos del CaseManifest:** candidate_rois, cohorte sin pedigree
individual (o pedigree poblacional si aplica).

**Falta por construir:** tool PLINK general (ya está el binario),
flujo de cohorte.

---

## Lo que comparten los cuatro

Todos empiezan con un `CaseManifest`. Todos usan el mismo motor de
grafo, las mismas tools básicas, el mismo mecanismo de trazabilidad.

## Lo que los diferencia

| Caso | pedigree | hpo | rois | bams | prioridad |
|------|----------|-----|------|------|-----------|
| 1    | ✓ (extraído) | ✓ (interrogado) | ✗→✓ (derivado) | ✗ | clínica |
| 2    | ✓ | ✓ | ✓ | ✓ | clínica |
| 3    | ✗ | ✗ | ✓ | ✓ | investigadora |
| 4    | ✗ o cohorte | ✗ | ✓ | ✓ | estadística |

## Prioridades de construcción

1. Carga de PDFs (caso 1).
2. HPO → genes con Phen2Gene (caso 1).
3. Cliente OMIM/Orphanet (caso 1).
4. Modo interrogación (capa agéntica, caso 1).
5. Procesamiento en lote (casos 3, 4).
6. Tool PLINK general (caso 4).

## Notas

- El modo interrogación debe hacer como máximo 3-5 preguntas
  agrupadas, ofrecer opciones cerradas cuando pueda, y permitir
  "exploración sin fenotipo" como salida de escape. Nunca insistir
  si el genetista dice "no tengo más información".
- Los cuatro casos cubren el espectro que Alkuraya conoce bien:
  clínica (1, 2), poblacional (3), estadístico (4).
