# De pipeline a decisión

Documento corto para fijar el encuadre del proyecto. No es diseño
técnico, es la distinción que define qué es PublicGenomicAgent y
qué no es.

## El pipeline genómico tradicional

    FASTQ -> BAM -> VCF -> VCF anotado -> (interpretación humana)

Todo el trabajo es transformación de datos. La decisión clínica vive
fuera del sistema, en la cabeza del genetista, y no queda registrada.
El resultado es un fichero, no una conclusión.

## PublicGenomicAgent

    BAM + consulta clínica -> decisión cuantificada

El trabajo empieza donde el pipeline termina. El BAM es el punto de
partida (el alineamiento es responsabilidad del laboratorio). La
consulta clínica es la entrada real. La decisión cuantificada es el
producto.

Componentes de la decisión:

- **Estadística mendeliana**: cada variante candidata tiene una
  probabilidad asociada basada en herencia, no en heurísticas.
  ¿De novo? ¿Recesivo? ¿Compuesto? ¿Consistente con el pedigree?
- **Pangenomización local**: corrección del reference bias para
  poblaciones infrarrepresentadas. Un VCF delta demuestra la ganancia.
- **Trazabilidad**: cada paso del proceso queda registrado, con
  hashes, versiones de binarios, y decisiones tomadas. El caso es
  auditable.

## Consecuencias de diseño

1. **El VCF es un medio, no un fin.** Es el lenguaje intermedio
   sobre el que se construye la decisión. Sin VCF no hay decisión,
   pero la decisión no es el VCF.

2. **La estadística mendeliana es un motor de primera clase.** No es
   un filtro auxiliar. Es lo que convierte observaciones en
   probabilidades. Su diseño merece la misma atención que la
   pangenomización local.

3. **El prompt clínico es entrada formal, no UI.** Cuando llegue la
   capa agéntica, el prompt no es una caja de chat. Es el input
   clínico: fenotipo, familia, sospecha, restricciones. Se parsea,
   se valida, y se registra.

4. **La cuantificación debe ser explícita.** Cada output del sistema
   debe poder responder: "¿con qué confianza?". No solo "aquí está
   la variante", sino "esta variante X tiene estas evidencias, con
   estos pesos, y estas son las alternativas consideradas".

5. **El caso es la unidad de trabajo, no el fichero.** Un caso
   incluye contexto, proceso, decisión y trazabilidad. No es un VCF
   suelto.

## Encuadre para comunicación externa

Frase corta:

> Yo no voy de FASTQ a VCF. Voy de BAM y una consulta clínica a
> una decisión cuantificada, con estadística mendeliana y
> pangenomización local, y todo el proceso queda auditable.

Esto es lo que diferencia el proyecto de un pipeline más. Es el
hueco real en genómica clínica hoy: no falta capacidad de procesar
datos, falta capacidad de producir decisiones reproducibles,
auditables y transferibles entre centros.
