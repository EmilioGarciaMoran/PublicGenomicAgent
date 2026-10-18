# PublicGenomicAgent — pitch

**A privacy-first, laptop-scale copilot for clinical genomics.**

Dado un caso clínico (fenotipo, familia, BAMs ya alineados),
PublicGenomicAgent completa el análisis genómico end-to-end en un
portátil de consumo. Sin clúster, sin GPU. Y con el LLM local por
defecto, sin red: el usuario puede abrir un proveedor externo bajo
consentimiento explícito, pero el caso de uso principal no lo
necesita. Cuando el pipeline lineal falla —duplicaciones
segmentarias— activa pangenoma local dirigido para recuperar las
variantes patogénicas que se pierden.

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

### 1. Análisis genético completo en un portátil de consumo

**No alineamos.** Partimos de BAMs ya alineados por el pipeline
habitual del hospital. El valor del sistema está en la capa de
interpretación dirigida por fenotipo, no en el alineamiento.

La entrada puede ser un **relato clínico** (que el sistema
estructura en un `CaseManifest` con validación estricta) o un
`CaseManifest` ya construido. En ambos casos, el sistema **no
ejecuta ninguna tool** hasta que el manifest está completo y
correcto (ver sección "Principio de diseño: manifest-first").

Dado un caso clínico, PublicGenomicAgent completa el análisis en
un ordenador de consumo:
