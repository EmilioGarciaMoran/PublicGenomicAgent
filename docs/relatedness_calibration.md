# Calibración de la verificación familiar con GIAB

La verificación familiar (`verify_relatedness`) no bloquea nunca un trío
verdadero y no da nunca por bueno a un padre falso, medido sobre exomas
reales de Genome in a Bottle (GIAB). Con menos de unos 200 SNV comunes
informativos, el resultado es «no verificable», no «aprobado».

Ejecutado el 2026-10-31 con `scripts/calibrate_relatedness_giab.py`,
`samtools` 1.19.2 y `bcftools` 1.19 (htslib 1.19.1).

## Datos

| Individuo | Papel en la calibración | Exoma |
|---|---|---|
| HG002 (NA24385) | hijo | Oslo University Hospital, Illumina, Agilent SureSelect v5, bwa-mem, GRCh37 |
| HG003 (NA24149) | padre | ídem |
| HG004 (NA24143) | madre | ídem |
| HG005 (NA24631, trío chino) | padre falso (otra familia) | ídem |

URLs y MD5 tomados de
[genome-in-a-bottle/giab_data_indexes](https://github.com/genome-in-a-bottle/giab_data_indexes).
No se descargaron los BAM completos: htslib leyó por HTTPS, a través del
índice, solo los cromosomas 21 y 22 (280–330 MB por muestra).

Panel de sitios: 3.471 SNV bialélicos de 1000 Genomes fase 3 con
frecuencia alélica entre 0,2 y 0,8, dentro de las 21.366 regiones
cubiertas a 20× o más en HG002. Los individuos GIAB no forman parte de
1000 Genomes, así que no hay circularidad.

## Parejas con todos los sitios

| Pareja | Relación real | n | IBS0 | IBS0/n | relatedness |
|---|---|---|---|---|---|
| HG003–HG002 | padre–hijo | 3.448 | 1 | 0,0003 | 0,553 |
| HG004–HG002 | madre–hijo | 3.459 | 7 | 0,0020 | 0,528 |
| HG003–HG004 | matrimonio, no emparentados (askenazíes) | 3.439 | 300 | 0,0872 | 0,083 |
| HG002–HG005 | no emparentados, otra ascendencia | 3.424 | 390 | 0,1139 | −0,197 |
| HG003–HG005 | no emparentados, otra ascendencia | 3.415 | 409 | 0,1198 | −0,189 |
| HG004–HG005 | no emparentados, otra ascendencia | 3.417 | 408 | 0,1194 | −0,200 |

- En padre/madre–hijo, IBS0 solo aparece por error de genotipado:
  entre un 0,03 % y un 0,2 % de los sitios.
- Entre no emparentados de la misma población endogámica (askenazí),
  IBS0 aparece en el 8,7 % de los sitios y el relatedness de fondo es
  0,08. Es la referencia más parecida a una población saudí que ofrece
  GIAB.
- Entre individuos de ascendencias distintas el relatedness sale
  negativo, como es habitual con este estimador. No afecta a las
  comprobaciones.

## Regla de decisión resultante

Para cada pareja declarada padre/madre–hijo:

- **Umbral de IBS0:** el menor recuento con probabilidad menor de 1e-3 si
  la pareja fuera realmente padre/madre–hijo, suponiendo un error de
  genotipado del 0,5 % (2,5 veces el peor observado). Superarlo es un
  fallo.
- **Potencia:** si con los sitios disponibles la probabilidad de detectar a
  un no-progenitor de la misma población (tasa de IBS0 del 8,7 %) es
  menor del 99,9 %, la pareja es «no verificable». En la práctica hacen
  falta unos 200 sitios informativos.

La regla fija anterior («falla si IBS0/n > 1 %») daba, con estos mismos
datos, un 1,2 % de falsos fallos con 100 sitios y un 4,4 % con 150: con
paneles pequeños bastaban dos errores de genotipado para cruzar el
umbral.

## Resultados por tamaño de panel

500 paneles simulados por tamaño, submuestreando los sitios informativos
(semilla 1).

**Trío verdadero** (debe pasar):

| Sitios del panel | passed | warning | failed | not_verifiable |
|---|---|---|---|---|
| 50 | 0 | 0 | 0 | 1,000 |
| 100 | 0 | 0 | 0 | 1,000 |
| 150 | 0 | 0 | 0 | 1,000 |
| 200 | 0,986 | 0,014 | 0 | 0 |
| 300 | 0,994 | 0,006 | 0 | 0 |
| 500 | 1,000 | 0 | 0 | 0 |
| 1.000 | 1,000 | 0 | 0 | 0 |

**Padre falso: HG005 declarado padre de HG002** (debe fallar):

| Sitios del panel | passed | warning | failed | not_verifiable |
|---|---|---|---|---|
| 50 | 0 | 0 | 0 | 1,000 |
| 100 | 0 | 0 | 0,994 | 0,006 |
| 150 | 0 | 0 | 1,000 | 0 |
| 200 | 0 | 0 | 1,000 | 0 |
| 300 | 0 | 0 | 1,000 | 0 |
| 500 | 0 | 0 | 1,000 | 0 |
| 1.000 | 0 | 0 | 1,000 | 0 |

Los avisos del trío verdadero (1,4 % con 200 sitios) proceden del
matrimonio HG003–HG004: con pocos sitios, el ruido lleva a veces su
relatedness de fondo (0,08) por encima de 0,3. Un aviso no bloquea el
análisis.

## Consecuencia para paneles

Un panel con menos de unos 200 SNV comunes cubiertos no permite
verificar la paternidad o maternidad con garantías. El resultado será
«no verificable» y solo se podrá continuar con una excepción humana
registrada (`--allow-unverified-relatedness`). Opciones para esos
paneles: aprovechar las lecturas fuera de diana o añadir unos cientos
de SNV de control al diseño.

## Límites

- El padre falso probado es de otra ascendencia. Para un no-progenitor
  de la misma población, la potencia está calculada con la tasa
  askenazí (8,7 %), no medida directamente.
- No se ha probado una madre falsa: GIAB no tiene una mujer no
  emparentada con exoma de Oslo.
- Solo se han usado los cromosomas 21 y 22 de un exoma de captura, en
  GRCh37. Los paneles reales pueden tener otra distribución de cobertura
  y de frecuencias.
- Las frecuencias del panel son globales de 1000 Genomes. Para cohortes
  saudíes conviene un panel con frecuencias de Oriente Medio.

## Cómo reproducirlo

Requiere red y el entorno `pga-hts`:

```bash
pga-env
python scripts/calibrate_relatedness_giab.py --workdir ~/pga_giab_calibration
```

Cada paso guarda su resultado en `--workdir` y se reaprovecha al
relanzar. El informe completo queda en `report.md` y `report.json`.
