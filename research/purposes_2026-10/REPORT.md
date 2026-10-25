# Propósitos publicables para relatipy — informe de la sesión del 2026-10-25

Este directorio no forma parte del paquete (los paquetes se construyen solo desde
`src/` y `MANIFEST.in` no lo incluye). Contiene los scripts, los resultados y las
notas de literatura que respaldan cada afirmación. Ningún archivo de la librería
fue modificado.

<!-- SECTION:PURPOSES -->

## 1. Qué hace la librería hoy (verificado leyendo el código y ejecutándolo)

- Geodésicas **timelike** en Kerr, Boyer–Lindquist, núcleo en C. Integradores
  `radau` (por defecto), `dop853`, `dp45` y `projection_radau`. La proyección
  fija la posición y corrige la 4-velocidad hacia la norma, E, Lz y Q iniciales.
- `KerrMcmcModel`: offsets en el cielo y v_los con retardo de Rømer resuelto de
  forma autoconsistente (`t_obs`). **No** propaga fotones (sin lente ni Shapiro).
  La posición en el cielo es la proyección de las coordenadas cartesianas BL
  (`x = sqrt(r²+a²) sinθ cosφ`, …). El redshift es `1+z = u^t + u_away`
  (`native/src/geodesic/physic/kerr_observables.c`, línea 143), es decir, un
  momento del fotón de espacio plano.
- El núcleo C ya contiene `rp_kerr_null_geodesic_rhs`, usado con argumento 1.0
  para el caso timelike (`native/src/geodesic/kerr_mcmc.c`). Exponer fotones
  está, por tanto, cerca.
- No admite fuerzas externas, masas extendidas ni eventos definidos por el
  usuario, y no libera el GIL (se paraleliza con procesos).

## 2. Resultados numéricos de esta noche

Sgr A*: M = 4.297e6 Msun, R0 = 8.277 kpc, 1 r_g = 5.124 μas = 0.0424 AU.
Todos los experimentos son sintéticos (inyección/recuperación). Los scripts
están en `scripts/` y sus salidas en `results/`.

### 2.1 Velocidad (`bench.py`, `speed_kg.py`, CPU casi libre)

| tarea | tiempo |
|---|---|
| S2, 200 épocas, 34 años, `t_eval` (dop853) | 6.5 ms (≈150 evals/s/núcleo) |
| ídem con inversión de Rømer (`t_obs`) | 8.6 ms |
| ídem con rtol de 1e-4 a 1e-13 | sin cambio apreciable (domina el muestreo) |
| `radau` en vez de `dop853` | ≈2.2× más lento |
| KerrGeoPy (analítico, Python): StableOrbit + 200 muestras en tiempo de Mino | 15.9 ms |
| KerrGeoPy + inversión t(λ) en 200 tiempos coordenados | 42.6 ms (sin Rømer, proyección ni RV) |
| Kerr exacto, 1000 órbitas a 300 r_g (QPE) | 0.25–0.39 s |

Con 4 núcleos, un MCMC de 10⁶ evaluaciones de S2 tarda ≈30 min. El Kerr exacto
en C cuesta lo mismo o menos que un integrador 1.5PN en numba (0.9 ms/órbita
frente a 0.25–0.39 ms/órbita).

### 2.2 Validación del límite post-newtoniano (`validate_pn.py`)

El término espín-órbita 1.5PN `a = (2/r³)[3(n·S)(n×v) + (v×S)]` reproduce la
precesión nodal del geodésico Kerr exacto a 1.3e-4 relativo (a = 3000 r_g,
e = 0.3, a* = 0.9, 20 órbitas). Las otras cuatro variantes de signo y orden
fallan por factores 2, ½ y −½. La fórmula analítica de Lense–Thirring queda a
0.9 %, como se espera de las correcciones de orden superior.

### 2.3 Señal irreducible de frame-dragging (`scan2.py`)

Se inyecta Kerr con a* = 0.9 y se reajusta con a* = 0, dejando libres los 8
parámetros. La campaña dura 20 años, con 4 épocas/año y 30 épocas densas en cada
pericentro. Residuo astrométrico rms (máximo), en μas:

| a [AU] | e | r_p [r_g] | rms (máx) |
|---|---|---|---|
| 1035 | 0.90 | 2441 | 0.37 (2.1) |
| 1035 | 0.976 | 586 | 3.7 (18.7) |
| 600 | 0.976 | 340 | 17.3 (68) |
| 400 | 0.976 | 226 | 32.7 (146) |
| 250 | 0.95 | 295 | 20.7 (121) |

S2 (r_p ≈ 2825 r_g) deja 0.1–0.4 μas rms. El espín es invisible en S2 con
cualquier instrumento previsto. Las filas con r_p < 100 r_g del JSON no son
fiables: el ajuste sin multi-arranque cae en mínimos locales.

### 2.4 El cuello de botella no es la dinámica, es el mapeo al cielo (`bias.py`, `bias_gauge.py`)

La verdad es Kerr exacto (BL). El modelo es 1PN o 1.5PN en gauge armónico
(numba), con la misma proyección y el mismo Rømer, solo astrometría y todos los
parámetros de ruido libres.

| caso | rms / máx [μas] | f_SP ajustado |
|---|---|---|
| S2, Newtoniano | 38.0 / 218 | — |
| S2, 1PN armónico | 0.91 / 3.66 | — |
| S2, 1PN con r → r + 1 r_g (BL) | **0.010 / 0.037** | — |
| S2, 1PN × f_SP, armónico | 0.20 / 0.40 | **1.0424** |
| S2, 1PN × f_SP, BL | 0.004 / 0.034 | **1.00045** |
| r_p = 424 r_g, 1PN, BL | 0.012 / 0.064 | 1.0028 |

- La truncación 2PN de la dinámica es despreciable (≤0.06 μas) incluso a
  r_p ≈ 400 r_g.
- Elegir coordenadas (armónicas o BL) para proyectar al cielo cambia la
  predicción en O(r_g) ≈ 5 μas, y ese cambio es degenerado con f_SP: sesga f_SP
  un 4.2 % en S2. GRAVITY 2024 da σ(f_SP) ≈ 0.11–0.13, así que hoy no importa,
  pero sí con la precisión de GRAVITY+ y el ELT.
- La astrometría relativista estándar (IAU/Gaia) se formula en gauge armónico:
  línea recta más deflexión O(r_g/b). Por eso es la proyección BL de
  `KerrMcmcModel` la que carga el término extra. Arreglo inmediato: transformar a
  coordenadas armónicas o isotrópicas antes de proyectar. Arreglo correcto:
  propagar fotones.
- Sesgo de espín de un modelo PN armónico (solo astrometría): S2 → a* = 2.39
  (verdad 0.9), sin sentido físico porque el ajuste usa el espín para absorber el
  gauge. Estrella con r_p = 424 r_g → 0.911. r_p = 302 r_g → 0.910.
- Los ajustes solo con astrometría deslizan M a lo largo de la degeneración
  M/D³ (verificado). Un análisis real debe incluir la RV.

### 2.5 Espín: Fisher y MCMC con Kerr exacto (`fisher*.py`, `mcmc_spin*.py`)

Precisión σ(a*) para a* = 0.9 con eje (θ, φ) = (60°, 30°), campaña de 20 años:

| estrellas | 50 μas / 10 km/s | 30 / 10 | 10 / 5 | 3 / 2 |
|---|---|---|---|---|
| S2 | 6.9 | 4.4 | 1.5 | 0.44 |
| S62-like (a = 749 AU, e = 0.976, r_p = 424 r_g) | 0.17 | 0.10 | 0.035 | 0.011 |
| S4716-like | 1.7 | 1.6 | 0.73 | 0.26 |
| S2 + S62 + S4716 + S29 | 0.126 | 0.078 | 0.027 | 0.008 |

Promediado sobre 60 ejes isótropos (S62-like sola): con 30 μas, σ(a*) mediana
0.27 (p16–p84: 0.10–0.48) y solo ≈15 % de las orientaciones bajo 0.1; con
10 μas / 5 km/s, mediana 0.09 (0.03–0.16), eje a 11–33°. El caso (60°, 30°) era
favorable.

MCMC (emcee, 17 parámetros, S2 + S62-like, 10 μas / 5 km/s). Con 40 walkers
(≈2.3 por dimensión) **el ensemble colapsó**: los anchos quedaron en 0.6× la cota
de Cramér–Rao y estables, lo que parecía convergencia. Con 96 walkers
inicializados desde la covarianza de Fisher se obtuvo a* = 0.949 ± 0.025,
cosθ = 0.35 ± 0.12 y φ = 18 ± 7°, frente a Fisher 0.035 / 0.137 / 9.8°: es
consistente, y el borde a* ≤ 1 explica el ancho algo menor. Lección: usar ≥4–6
walkers por dimensión e inicializar desde Fisher.

**S301** es real: Nature 2026, arXiv:2607.12664, ESO eso2612. K ≈ 19.3,
P ≈ 8.7 años, e ≈ 0.982, r_p ≈ 280 r_g, próximo pericentro en 2031. El Fisher
siguiente usa a* = 0.7, 60 ejes isótropos, los ángulos orbitales de S2 como
sustitutos de los de S301 y 4 épocas/año más 30 por pericentro:

| campaña | σ_ast / σ_RV | σ(a*) mediana (p16–p84) | eje |
|---|---|---|---|
| 2026–2035 (1 pericentro) | 100 μas / 1 km/s (supuestos de Nature) | 0.60 (0.27–1.02) | 91° |
| 2026–2035 | 10 μas / 1 km/s | 0.21 | 32° |
| 2026–2040 (2 pericentros) | 100 μas / 1 km/s | 0.19 (0.12–0.28) | 24° |
| 2026–2040 | 30 μas / 10 km/s | 0.33 | 46° |
| 2026–2040 | 10 μas / 1 km/s | 0.079 (0.04–0.13) | 11° |

Con un solo pericentro no basta. La RV a 1 km/s aporta la mayor parte de la
información.

### 2.6 Redshift a precisión MICADO (`rv_term.py`, `rv_s301.py`)

En el pericentro de S301 (|v_los| máx ≈ 23 000 km/s), diferencia bruta, antes de
reajustar, respecto de la referencia Doppler de marco estático
`1+z = u^t + ĝ·u_ON/√(1−2/r)` sin deflexión:

| modelo | error máx |
|---|---|
| 1PN aditivo (v_los + v²/2 + M/r) | 256 km/s |
| 1PN + término cruzado | 100 km/s |
| SR × GR multiplicativo, fotón plano | 84 km/s |
| `KerrMcmcModel` actual | 77 km/s |

El término métrico que omite la librería vale 2.4 km/s en S2, 42 km/s a
r_p = 424 r_g y 70 km/s a 300 r_g. La deflexión de la luz aporta otro
O(r_g/r · v) ≈ 80 km/s, así que la propia referencia es incierta en esa
cantidad. Conclusión: con RV a 1 km/s, S301 exige trazar geodésicas nulas en Kerr
dentro del MCMC. Falta cuantificar cuánto absorbe un reajuste.

### 2.7 QPE (`qpe.py`, `qpe2.py`)

Un orbitador cruza un disco ecuatorial. M = 4e5 Msun. Se comparan los tiempos de
cruce con Rømer: verdad Kerr exacto; modelo 1.5PN en gauge BL con todos los
parámetros y a* libres.

| órbita | ventana del ajuste | rms / máx | a* ajustado (verdad 0.7) |
|---|---|---|---|
| a = 300 r_g, e = 0.2 (P = 17.9 h) | 100 órbitas | 0.088 s / 0.43 s | 0.652 |
| a = 300 r_g | 300 órbitas | 0.024 s / 0.09 s | 0.669 |

<!-- SECTION:QPE60 -->

A 300 r_g el 1.5PN basta para el timing: los tiempos de erupción se miden con
errores de ~10² s. El espín sí sale sesgado (−0.03 a −0.05).

### 2.8 Integradores y la proyección (`integ_bench2.py`)

Órbita a* = 0.9, p = 10, e = 0.5, x = 0.5, 1000 órbitas radiales, pasos
aceptados, comparación con la solución analítica de KerrGeoPy:

| método | rtol | tiempo | dE/E | dLz/Lz | dQ/Q | \|Δφ\|máx |
|---|---|---|---|---|---|---|
| dop853 | 1e-12 | 0.53 s | 3.5e-12 | 2.0e-9 | 9.3e-10 | 1.6e-5 |
| radau (por defecto) | 1e-12 | 6.3 s | 9.1e-10 | 2.0e-7 | 2.2e-7 | 5.8e-3 |
| projection_radau | 1e-12 | 6.2 s | 4.6e-16 | 4.9e-16 | 8.9e-16 | **1.9e-6** |
| projection_radau | 1e-10 | — | falla en τ = 1.47e5 de 3.48e5 ("invariant projection failed; native status 10") | | | |
| dp45 | 1e-12 | 1.9 s | 2.5e-10 | 4.6e-9 | 7.1e-9 | 6.5e-4 |

- La proyección lleva los cuatro invariantes a precisión de máquina y mejora la
  fase de radau ≈3000×.
- `radau` es disipativo (L-estable, no simplético): es ~12× más lento y ~400×
  peor en fase que dop853 al mismo rtol. Es mala opción por defecto para
  geodésicas no rígidas.
- `projection_radau` falla con rtol de 1e-8, 1e-9 y 1e-10, siempre en
  r − r_p = +0.006…0.013 r_g, es decir, justo en el punto de retorno radial. Con
  la posición fija y el error acumulado, R(r) < 0 y no existe una 4-velocidad
  real que cumpla las cuatro restricciones.
- Las derivas de E, Lz y Q evaluadas en muestras `tau_eval` están dominadas por
  la interpolación (~1e-5), no por el integrador.

### 2.9 Robustez (`robust.py`)

150 perturbaciones aleatorias (10⁻³) por caso. Con `t_obs` hubo 0 fallos en todos
los casos. Con `t_eval` + dop853 hubo 1/150 (S4716-like) y 3/150 (S29); con
radau, 0. Además, el Fisher encontró un fallo de `t_obs` con rtol de 1e-13 y
1e-12 para S4716-like con M·(1 − 1e-5). En un MCMC hay que devolver −∞ ante
`IntegrationError` y medir la tasa de fallos.

<!-- SECTION:BUGS -->

## 5. Bloqueos de esta sesión

La política de red del entorno denegó arxiv.org, export.arxiv.org,
cdsarc.cds.unistra.fr, vizier.cds.unistra.fr, ui.adsabs.harvard.edu,
iopscience.iop.org, aanda.org, doi.org y zenodo.org. Por eso:

- la literatura se basa en resúmenes de búsqueda; cada cita está marcada en
  `literature/` como verificada o solo snippet;
- no se descargaron los datos públicos de S2 (Gillessen+2017, CDS J/ApJ/837/30,
  `table5.dat`; Chu+2018 J/ApJ/854/12).

Para arreglarlo: menú del entorno cloud (barra de título de la sesión) → Edit →
Network access → añadir esos dominios en *Allowed domains*, dejando marcada la
casilla de gestores de paquetes. Documentación:
https://code.claude.com/docs/en/cloud-environments#network-access

## 6. Reproducir

```console
python -m pip install . emcee kerrgeopy numba      # desde la raíz del repo
cd research/purposes_2026-10/scripts
python validate_pn.py                 # §2.2
python scan2.py                       # §2.3 (≈11 min)
python bias.py S2 S62-like            # §2.4
python bias_gauge.py gauge S2 S62-like
python fisher.py 0.9; python fisher_axes.py; python fisher_s301.py   # §2.5
python mcmc_spin2.py 4000             # §2.5 (≈32 min en 4 núcleos)
python rv_s301.py                     # §2.6
python qpe.py 300 0.2 40 0.7 1000 100,300    # §2.7
python integ_bench2.py 0.9 10 0.5 0.5 1000   # §2.8
python robust.py 150                  # §2.9
```

<!-- SECTION:REFS -->
