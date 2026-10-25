# Propósitos publicables para relatipy — informe de la sesión del 2026-10-25

Este directorio no forma parte del paquete (los paquetes se construyen solo desde
`src/` y `MANIFEST.in` no lo incluye). Contiene los scripts, los resultados y las
notas de literatura que respaldan cada afirmación. Ningún archivo de la librería
fue modificado.

## 0. Veredicto

- **S2:** su espín no se puede medir con ningún instrumento previsto. La señal
  irreducible de frame-dragging es de 0.1–0.4 μas rms y σ(a*) ≥ 1.5 incluso con
  10 μas. Además, los sistemáticos de propagación de luz (4–25 μas) son mayores
  que la señal.
- **S301** (Nature 2026; r_p ≈ 280 r_g; pericentros en 2031 y ~2040) es el
  objetivo correcto. Ahí la dinámica exacta importa (2PN ≈ 44 μas/órbita) y el
  modelo de RV actual se equivoca entre 3 y 600 km/s cuando MICADO apunta a
  1 km/s.
- La biblioteca es útil hoy como **referencia exacta y rápida**. Para publicar
  ciencia con S301 le falta propagar la luz, de forma barata con una función T
  de primer orden validada con rayos exactos.

## 3. Los cinco propósitos, ordenados por publicabilidad × necesidad real de Kerr exacto y rápido

### P1. Presupuesto de sistemáticos del modelo observable, de S2 a S301
- **Problema real.** GRAVITY mide f_SP, el redshift f y masas extendidas de
  <1200 Msun con modelos PN que proyectan coordenadas al cielo. A precisión de
  GRAVITY+, ELT o MICADO, ¿cuánto sesgo introduce la propagación de luz
  omitida?
- **Evidencia de esta noche.**
  - f_SP se desplaza un 4.2 % en S2 entre mapas de coordenadas.
  - La recta armónica erra 2M·tan(ψ/2), entre 4 y 25 μas en S2.
  - La recta BL erra M·sinψ·tan²(ψ/2).
  - El error de RV es (M/r)·v_los·tan²(ψ/2), de 3 a 609 km/s en S301.
  - La dinámica 2PN es despreciable en S2 pero no en S301.
- **Por qué hace falta la biblioteca.** Para generar mocks con Kerr exacto y
  rayos trazados, y reajustarlos miles de veces con un modelo tipo GRAVITY
  (velocidad).
- **Falta.** Exponer el RHS nulo del núcleo C; implementar el observable T de
  primer orden; mocks trazados con rayos, reajustados, que den Δf_SP, ΔR0 y
  Δa* en unidades de σ.
- **Riesgo de novedad:** moderado (Grould+2017; Zhang, Lu & Yu 2015). **Revista:**
  A&A o MNRAS.

### P2. S301: forecast independiente del vector de espín de Sgr A*
- **Problema real.** Nature 2026 anuncia el espín "en ~una década" con 100 μas
  y 1 km/s.
- **Evidencia.** Fisher con ángulos sustitutos: con un pericentro, σ(a*)
  mediana 0.60; con dos, 0.19; con 10 μas / 1 km/s, 0.079. La RV a 1 km/s
  aporta casi toda la información, y es justo donde el modelo actual falla. El
  MCMC con Kerr exacto valida el Fisher (S2 + S62-like).
- **Falta.**
  - Elementos publicados de S301.
  - Espín cartesiano.
  - Parámetros de marco (x0, y0, vx0, vy0, vz0).
  - Masa extendida.
  - Datos de Asimov o muchas realizaciones.
  - El observable de P1.
  - Reproducir antes exactamente el setup del paper Nature.
- **Riesgo:** alto y sensible al tiempo (el equipo de Nature/GRAVITY, Zhang/Yu,
  Iorio). Lo diferencia el código abierto con dinámica exacta, el modelo de luz
  validado y el presupuesto de sesgos. **Revista:** A&A o ApJ.

### P3. S301 en un cúmulo realista: espín frente a masa extendida y perturbadores
- **Problema real.** Merritt+2010 sostiene que las perturbaciones del cúmulo
  enmascaran el frame-dragging más allá de ~0.5 mpc; S301 tiene a ≈ 3.3 mpc. La
  masa extendida (<1200 Msun, GRAVITY 2024) es degenerada con el término de
  espín en el plano para una sola estrella.
- **Falta (extensión de la biblioteca).** Fuerzas externas: potencial de masa
  extendida y N perturbadores (Newtonianos sobre el fondo de Kerr). Después,
  mocks con perturbadores reajustados por espín.
- **Riesgo:** moderado; es Merritt+2010 rehecho con Kerr exacto, la órbita real
  y los límites actuales. **Revista:** MNRAS o PRD.

### P4. Verosimilitud diferenciable para geodésicas de Kerr (métodos)
- **Problema real.** Los ajustes multi-estrella con perturbadores tienen más de
  30 parámetros. Con emcee el ensemble colapsa, como pasó aquí con 40 walkers.
- **Propuesta.** Ecuaciones variacionales en C que den gradientes exactos de
  los observables, para HMC/NUTS y matrices de Fisher exactas. Incluir los
  integradores con proyección (invariantes a 1e-16, fase ×3000 mejor que
  radau), con el fallo en puntos de retorno resuelto.
- **Revista:** JOSS (software) más A&C, ApJS o CPC (método). **Riesgo:** bajo.
  El revisor tropezará con los bugs de la §4: hay que arreglarlos antes.

### P5. QPE de órbita cercana: referencia exacta para los modelos de timing
- **Problema real.** El timing de QPE (GSN 069, eRO-QPE2…) usa geodésicas
  expandidas a 3PN (QPE-FIT, Chakraborty+2025) o esquemas osculantes
  (Zhou+2024).
- **Evidencia.** A 300 r_g el 1.5PN basta: 0.02–0.09 s en 300 órbitas. A
  60 r_g falla: **225 s rms y 544 s de máximo en 300 órbitas, con a* = −0.15**
  frente a 0.7.
- **Falta.** Sesgo/σ para a = 10–100 r_g con jitter realista (10²–10³ s), un
  disco desalineado que precesa por LT y Rømer; comparar con el 3PN de QPE-FIT.
- **Riesgo:** muy alto (Zhou+ podría usar ya Kerr exacto; sin verificar). Las
  soluciones analíticas en tiempo de Mino son más rápidas para geodésicas puras.

### Descartado
Hot spots/flares a 6–12 r_g: el campo está muy poblado (GRAVITY/GYOTO,
Yfantis+2024, Aimar+2023) y la biblioteca no tiene trazado de rayos, emisión ni
polarización.

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

Con 4 núcleos, un MCMC de 10⁶ evaluaciones de S2 tarda ≈30 min. Comprobé que
`rtol` sí llega al núcleo C: las salidas cambian 0.05 μas con 1e-10, 4 μas con
1e-7 y 1.5 mas con 1e-4, y los pasos pasan de 4021 a 779. El tiempo plano se
debe al coste por época. KerrGeoPy es Python puro, así que la comparación de
velocidad no es justa frente a un código analítico compilado. El Kerr exacto
en C cuesta lo mismo o menos que un integrador 1.5PN en numba (0.9 ms/órbita
frente a 0.25–0.39 ms/órbita).

### 2.2 Validación del límite post-newtoniano (`validate_pn.py`)

El término espín-órbita 1.5PN `a = (2/r³)[3(n·S)(n×v) + (v×S)]` reproduce la
precesión nodal del geodésico Kerr exacto a 1.3e-4 relativo (a = 3000 r_g,
e = 0.3, a* = 0.9, 20 órbitas). Las otras cuatro variantes de signo y orden
fallan por factores 2, ½ y −½: el signo y el coeficiente quedan validados. La
coincidencia a 1e-4 es en parte casual. El cuadrupolo de Kerr aporta ≈ −0.7 % a
esta inclinación, (3/4)·a*·cos i/√p, y no está en el modelo PN. Para
separarlo hay que correr con +a* y −a*.

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

- **Corrección tras la revisión adversarial.** La primera lectura (que la
  proyección BL era la sesgada) estaba al revés. Lanzando rayos exactos en
  Schwarzschild (`scripts/chk.py`) se obtiene imagen = x_armónico,⊥ +
  2M·tan(ψ/2)·b̂ (ψ = ángulo observador–BH–estrella; error < 0.02 M hasta
  ψ = 120°). Errores frente a la verdad para ψ = 30/60/90/120/150°:
  - recta armónica: 0.54 / 1.16 / 2.0 / 3.45 / 7.2 M;
  - recta BL (relatipy): 0.04 / 0.29 / 1.0 / 2.6 / 6.7 M.

  La deflexión es O(M/b) pero desplaza la imagen O(M), el mismo orden que el
  cambio de coordenadas. Solo el mapa de coordenadas junto con la propagación
  de luz es invariante. Por tanto, el 4.2 % en f_SP mide la **propagación de luz
  omitida**: un modelo de recta armónica erraría 4–25 μas en S2, cuyo ψ recorre
  44–136°. GRAVITY 2024 da σ(f_SP) ≈ 0.11–0.13, así que hoy no domina; con
  GRAVITY+ y el ELT sí.
- La dinámica 2PN es despreciable **solo para S2**: la integral de Darwin da
  3.1e-6 rad/órbita, unos 0.7 μas. Para **S301** son **2.7e-4 rad/órbita, 44 μas
  en apocentro**, ≈14 % del término de Lense–Thirring en el plano para a* = 0.7.
  Mis ajustes solo con astrometría lo absorbían deslizando M ∝ R0³; con RV
  sesgará el espín. Para S301 la dinámica exacta sí importa.
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

**Corrección tras la revisión adversarial.** Mi primera referencia (Doppler de
marco estático sin deflexión) **no es una verdad válida**: le falta la
dirección real del fotón, y su error es −(2M/r)·tan(ψ/2)·(u·b̂). El ranking
"1PN 256 > … > relatipy 77 km/s" queda retirado. Rayos exactos en el pericentro
de S301 (r = 280 M, v ≈ 25 000 km/s), `scripts/chk.py`:

| ψ | 30° | 60° | 90° | 120° | 150° |
|---|---|---|---|---|---|
| relatipy − exacto [km/s] | 3 | 26 | 90 | 233 | 609 |
| referencia estática − exacto [km/s] | −42 | −53 | −1 | +155 | +563 |

En forma cerrada, el error de RV de relatipy en el pericentro es
(M/r)·v_los·tan²(ψ/2), con el mismo factor tan² que en la astrometría. MICADO
apunta a 1 km/s. Para ψ ≤ 120° basta una función de transferencia temporal de
primer orden T (Rømer + Shapiro, con 1+z = u^t + u^i ∂_i T y la imagen desde
∇⊥T), que da ≤0.02 M y ≤1.5 km/s. Los rayos nulos exactos sirven para validar
y para el lado lejano; no hacen falta dentro del MCMC. Las comparaciones sin
reajuste exageran el sesgo: el siguiente paso es reajustar o usar el formalismo
de Cutler–Vallisneri. Las correcciones de Kerr (lapso, g_tφ ~0.5 km/s;
deflexión por espín ≲1 km/s) cuentan para un forecast a 1 km/s.

### 2.7 QPE (`qpe.py`, `qpe2.py`)

Un orbitador cruza un disco ecuatorial. M = 4e5 Msun. Se comparan los tiempos de
cruce con Rømer: verdad Kerr exacto; modelo 1.5PN en gauge BL con todos los
parámetros y a* libres.

| órbita | ventana del ajuste | rms / máx | a* ajustado (verdad 0.7) |
|---|---|---|---|
| a = 300 r_g, e = 0.2 (P = 17.9 h) | 100 órbitas | 0.088 s / 0.43 s | 0.652 |
| a = 300 r_g | 300 órbitas | 0.024 s / 0.09 s | 0.669 |

| a = 60 r_g, e = 0.1 (P = 1.6 h) | 100 órbitas | 19.6 s / 36.8 s | 0.38 |
| a = 60 r_g | 300 órbitas (~20 días) | **225 s / 544 s** | **−0.15** |

A 300 r_g el 1.5PN basta para el timing: los tiempos de erupción se miden con
errores de ~10² s. El espín sí sale sesgado (−0.03 a −0.05). A 60 r_g el 1.5PN
falla por completo: los residuos son de minutos y el espín sale sin sentido.
Queda abierto si el 3PN de QPE-FIT basta; relatipy sirve de referencia exacta
para medirlo.

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

## 4. Bugs y limitaciones de la biblioteca, por importancia

1. **Redshift sin factores métricos ni momento real del fotón**
   (`kerr_observables.c:143`). Error (M/r)·v_los·tan²(ψ/2): 3–609 km/s en el
   pericentro de S301 y ~2 km/s en S2. Usar 1+z = −k·u/E, con k obtenido de T o
   de rayos.
2. **Sin propagación de luz**: ni Shapiro ni deflexión; el cielo es una
   proyección recta de BL. Error de imagen M·sinψ·tan²(ψ/2), hasta ~7 M
   (≈35 μas) en el lado lejano. Arreglo: T de primer orden con r_H = r − M.
3. **La inversión de Rømer omite Shapiro**: 100–200 s de retardo variable,
   ≈2 μas en S301.
4. **`radau` por defecto**: a igual rtol es ~12× más lento y ~400× peor en
   fase que dop853, y es disipativo. Conviene dop853 por defecto, con
   proyección opcional.
5. **`projection_radau` falla en el punto de retorno radial** con rtol ≥ 1e-10:
   el jacobiano es singular en u^r = 0 y no hay solución cuando R(r) < 0.
   Proyectar en el espacio de fases completo o usar variables de ángulo.
6. **Fallos esporádicos de `IntegrationError`** en `KerrMcmcModel`: ~1–2 % con
   `t_eval` + dop853 en órbitas largas excéntricas; raros con `t_obs` a rtol
   1e-12/1e-13.
7. **dLz deriva ~600× más que dE** con dop853: revisar la formulación y el
   redondeo.
8. **Espín**: falta una opción cartesiana y un manejo correcto del borde
   a* ≤ 1.
9. **`KerrMcmcModel` sin parámetros de marco**, masa extendida, factores de
   test (f_SP) ni fuerzas externas.
10. **RHS nulo no expuesto.**

`rtol` sí se transmite correctamente: verificado, no es bug.

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

## 7. Referencias clave

Estado: "s" significa que solo se vio en resúmenes de búsqueda porque arXiv y
ADS estaban bloqueados; "m" significa que la aportó el revisor de memoria.
**Verificar todas antes de citar.** Los detalles por cita están en
`literature/lit_*.md`.

- S301: Nature 2026, DOI 10.1038/s41586-026-10894-w, arXiv:2607.12664;
  ESO eso2612 (s).
- GRAVITY 2018 A&A 615 L15 (arXiv:1807.09409); 2020 A&A 636 L5
  (arXiv:2004.07187); 2022 A&A 657 L12 (arXiv:2112.07478); 2024 A&A 692 A242
  (arXiv:2409.12261) (s).
- Grould, Vincent, Paumard & Perrin 2017, A&A 608 A60 (arXiv:1709.04492) (s).
- Zhang, Lu & Yu 2015, ApJ 809 127: GR completa con fotones y MCMC de espín (m).
- Merritt, Alexander, Mikkola & Will 2010, PRD 81 062002 (arXiv:0911.4718) (s).
- Will 2008, ApJL 674 L25 (arXiv:0711.1677) (s).
- Waisberg et al. 2018 (arXiv:1802.08198) (s).
- QPE: Zhou et al. 2024 (arXiv:2401.11190, 2405.06429, 2411.18046); Chakraborty
  et al. 2025, QPE-FIT (arXiv:2508.20162); Arcodia et al. 2026
  (arXiv:2604.09788) (s).
- Software: KerrGeoPy (arXiv:2406.01413); PyGRO (arXiv:2504.20152);
  GC-OrbitFit (github pmplewa); GYOTO (Vincent+2011) (s).
- Flares: GRAVITY 2018 A&A 618 L10; Yfantis et al. 2024 (arXiv:2408.07120) (s).
