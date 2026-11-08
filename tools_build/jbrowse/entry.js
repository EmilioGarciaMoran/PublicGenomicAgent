// Visor JBrowse embebido para los informes de PGA.
// Recibe de Python (en <script id="pga-jbrowse-data">) la configuración
// completa de JBrowse y los ficheros en base64. Los ficheros se convierten
// en File y se registran como BlobLocation: JBrowse los lee de memoria,
// sin red y abriendo el HTML con file://.
import React from 'react'
import { createRoot } from 'react-dom/client'
import { createViewState, JBrowseLinearGenomeView } from '@jbrowse/react-linear-genome-view2'
import { setBlobMap } from '@jbrowse/core/util/tracks'

function b64ToFile(b64, name) {
  const bin = atob(b64)
  const bytes = new Uint8Array(bin.length)
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i)
  return new File([bytes], name)
}

const data = JSON.parse(document.getElementById('pga-jbrowse-data').textContent)
const blobs = {}
for (const [name, b64] of Object.entries(data.files)) blobs[name] = b64ToFile(b64, name)
setBlobMap(blobs)

const state = createViewState({ assembly: data.assembly, tracks: data.tracks, location: data.locus })
for (const id of data.show) state.session.view.showTrack(id)
window.__pgaJBrowse = state
createRoot(document.getElementById('pga-jbrowse')).render(
  React.createElement(JBrowseLinearGenomeView, { viewState: state }),
)
