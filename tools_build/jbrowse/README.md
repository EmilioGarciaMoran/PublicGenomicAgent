# Componente JBrowse embebido (E3)

El informe HTML incluye el visor lineal de JBrowse 2 compilado en un único
fichero: `src/publicgenomicagent/agent/vendor/jbrowse-lgv.min.js`, con su
hash en `vendor/SHA256SUMS` y las licencias de terceros en
`jbrowse-lgv.min.js.LEGAL.txt`. JBrowse es Apache-2.0; React es MIT.

No hace falta Node para usar PGA. Solo para recompilar:

    cd tools_build/jbrowse
    npm ci
    npm run build
    cd ../../src/publicgenomicagent/agent/vendor && sha256sum jbrowse-lgv.min.js > SHA256SUMS

`entry.js` lee el JSON de `#pga-jbrowse-data` (ensamblaje, pistas, locus y
ficheros en base64) y monta el visor en `#pga-jbrowse`, sin red.
