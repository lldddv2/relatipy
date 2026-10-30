# Notas de validación de relatipy (Kerr, partículas masivas)

Documento de traspaso para retomar en una sesión de código. Escrito el 2026-10-08.

Objetivo: una sección de validación al estilo de PyGRO (Della Monica 2025, §3), pero en Kerr
y solo con geodésicas temporales. Cada validación compara la integración de relatipy contra
un resultado teórico con referencia citada.

Leyenda de confianza usada abajo:

- **[num]** comprobado numéricamente hoy con relatipy 1.2.0 (código fuente local).
- **[fuente]** leído hoy directamente en la referencia.
- **[memoria]** citado de memoria; revisar en la referencia antes de publicarlo.

## 0. Estado

| # | Validación | Estado | Notebook |
|---|---|---|---|
| 1 | ISCO prógrado y retrógrado | **hecha** (2026-10-12): radio, frecuencia epicíclica, caída Mummery–Balbus; estabilidad de órbitas circulares vs PyGRO (2026-10-16, §1.4) | `01-…`, sección 1 |
| 2 | Arrastre de marco: caída ecuatorial desde el reposo en infinito | **hecha** (2026-10-12) | `01-…`, sección 2 |
| 3 | Precesión exacta (pericentro, ecuatorial) | **hecha**; barridos en $a$ y $e$ con $\chi=0.9$ | `01-…`, sección 3 |
| 4 | MCMC con S2 y S301 | por hacer; faltan datos de S301 | — |
| 5 | Kerr + disco delgado: campo (Foschi+ 2026), precesión apsidal coplanar (Tomaselli & Caputo + ⟨Φ⟩ de Piran+ 2026), nodal LT (Piran eq. 4) y disco (Piran eqs. 1, 10, 11, Tabla 1) | **hecha** (2026-10-25); fuentes en brain INF-037, revisión pendiente | `01-…`, sección 5 |

Sección 5 (2026-10-25):

- Posiciones 3D desde `PerturbedKerrMcmcModel.get_ra_dec_vr` con `t_eval` (sin Rømer) y dos observadores (permutación cíclica de ejes). Pasar `vec` también con χ=0: si no, el eje polar BL cae en el plano del disco y la integración falla.
- Requiere `relatipy._force_core` compilado. Se ejecutó con una copia limpia de `HEAD` (028d85b) en `PYTHONPATH`; el árbol tenía cambios en curso de otra sesión en `native/`.
- Sustitución propia (marcada en el notebook): promedio azimutal de $\Delta^{-1}$ = Foschi eq. 7 con $z^2\to z^2+H^2$ (Δ de Piran).
- Hallazgo: S301 con ω=0° se aparta 3–8 % de la eq. 1 de Piran; escala como $r_g/p$ (acoplamiento GR × disco cerca del pericentro).

Estado del notebook al 2026-10-12 (reemplaza lo que digan §2–§4 sobre implementación):

- Spin escrito $\chi$ en todo el notebook; $a$ es semieje mayor.
- Funciones compartidas: `equatorial_orbit`, `upward_crossings` (cruce con ajuste cúbico; reemplaza
  a `crossing_phases`), `glue_panels`.
- Integrador `projection_radau` en secciones 1 y 3. Sección 2 usa `radau`: `projection_radau`
  falla cerca de $r_+$ con $\chi\ge0.9$ (`invariant projection failed; native status 10`).
- Frecuencia epicíclica: perturbación $\delta=10^{-4}(1-r_{\rm ISCO}/r)$ y comparación contra
  $\Omega_r(\bar r)$ en el radio medio; error $10^{-7}$–$10^{-5}$. Con $\delta$ fijo domina el
  corrimiento por amplitud finita.
- Caída desde el ISCO: arrancar en $r_I(1-10^{-2})$; con $\epsilon=10^{-6}$ el término
  $\epsilon^3$ queda bajo el redondeo de $E,L_z$. Error en $\tau$ ~$10^{-7}$.
- `tau_eval` falla cerca del ISCO (reconstrucción interpolada, `native status 3`); usar `max_step`.
- Precesión: solo órbitas ecuatoriales; en inclinadas la modulación polar sesga la media ($1/N$).

Sección 1.4 (2026-10-16): análogo Kerr de PyGRO §3.2/Fig. 7. Radios $r_0/r_{\rm ISCO}=5/6,1,5/3$, $\chi=0.9$,
$\delta=10^{-4}$, 20 vueltas, `radau` `rtol=1e-12`, PyGRO `goal=12` `hmax`, KerrGeoPy `kg.Orbit` como referencia
analítica (r en los tiempos t de cada integrador). Con $\delta=10^{-8}$: KerrGeoPy resuelve la raíz doble solo a ~1e-8
(retrógrada estable sale "plunging") y la inestable retrógrada atraviesa la barrera por error de integración y cae. Fórmulas: BPT 1972 eqs. 2.9a, 2.10, criterio p. 351
(brain INF-031, revisión pendiente). Lyapunov no incluido: Pradhan 2012 eq. 104 tiene errata y Cardoso 2009 eq. 19
requiere combinar ecuaciones; decisión pendiente del usuario.

- **`projection_radau` desvía órbitas circulares retrógradas**: con $\delta=0$ y 5 vueltas, $|r/r_0-1|\sim4\times10^{-7}$
  (ISCO) y $\sim3\times10^{-8}$ ($5/3\,r_{\rm ISCO}$), igual con `rtol=1e-10` y `1e-13`; `radau`, `dop853` y PyGRO
  $\sim10^{-15}$. Prógradas no afectadas. Sin investigar en C; puede inflar errores retrógrados de §1.2–1.3.
  **Causa (2026-10-16)**: `rp_kerr_integrator_project` (`native/src/geodesic/integrators/kerr.c`) rehace
  $u^r=\pm\sqrt{\text{radial\_potential}}$ desde la normalización. Cerca de un punto de retorno (o en órbita circular)
  ese potencial es ~0 por cancelación de términos $O(1)$ con error ~$10^{-16}$; su raíz mete $|u^r|\sim10^{-8}$
  (medido: 1e-8 con proyección vs 1e-16 con `radau`). Independiente de `rtol`.
  **Corregido (2026-10-16)**: proyección ortogonal (nota `thesis/docs/v1/implementacion/2026-10-16-proyeccion-ortogonal.md`).
  Mismo tiempo (mediana 0.97×), error 2×–6000× menor; §1.4 ahora Relatipy ≈ PyGRO frente a KerrGeoPy. La caída
  ecuatorial de §2 con `projection_radau` ya llega a $r_+(1+10^{-3})$ para χ=0.5, 0.9, 0.998 (antes fallaba con χ≥0.9).
  §2 queda con `projection_radau` (decisión del usuario 2026-10-16), aunque ahí `radau` es ~7× más preciso
  ($\varepsilon_\phi$ ~5e-8 vs 7e-9); sigue ~10× mejor que PyGRO.

Cosas del entorno que hay que saber:

- `!pip install -Uq relatipy` (celda 0 del notebook) instala la **0.3.0** de PyPI. Todo lo de
  estas notas se corrió con el código local (1.2.0, `requires-python >= 3.11`). La API usada
  aquí (`Kerr.orbit(p=, e=, x=)`, `orbit.get_E()`, `sol.uR`, …) no se probó en 0.3.0.
- El notebook `01-…` ya fue reorganizado a mano después de añadir el ejemplo (8 celdas al momento
  de escribir esto); las funciones clave son `exact_precession` y `crossing_phases`.
- Hay un archivo vacío `.git/stale-index.lock.claude` que se puede borrar (un `git status`
  dejó un `index.lock` que no se pudo eliminar desde la sesión remota y se renombró).
- Este `.md` está dentro de `docs/`; si Sphinx/myst lo recoge, excluirlo en `conf.py` o moverlo.

## 1. Convenciones

- Unidades geométricas $G=c=M=1$; $a$ es el spin adimensional ($0\le a\le1$), coordenadas de
  Boyer–Lindquist, $\Delta=r^2-2r+a^2$, $\Sigma=r^2+a^2\cos^2\theta$, $r_\pm=1\pm\sqrt{1-a^2}$.
- Signo superior = órbita prógrada, inferior = retrógrada.
- Constantes de movimiento como las devuelve relatipy: `get_E()` $=-u_t$, `get_Lz()` $=u_\phi$,
  `get_Q()` = constante de Carter $Q$ (no $K$); todas con $M=\mu=1$.
- Potenciales (Carter 1968; forma de Fujita & Hikida 2009, §2) **[fuente]**:
  - $R(r)=[E(r^2+a^2)-aL_z]^2-\Delta\,[r^2+(L_z-aE)^2+Q]$, con $\Sigma\,\dot r=\pm\sqrt{R}$
  - $\Theta(\cos\theta)=Q-[Q+a^2(1-E^2)+L_z^2]\cos^2\theta+a^2(1-E^2)\cos^4\theta$
  - $\Sigma\,\dot\phi=\dfrac{a}{\Delta}[E(r^2+a^2)-aL_z]-aE+\dfrac{L_z}{\sin^2\theta}$
  - $\Sigma\,\dot t=\dfrac{r^2+a^2}{\Delta}[E(r^2+a^2)-aL_z]-a(aE\sin^2\theta-L_z)$
- **relatipy pide velocidades en tiempo coordenado** (`vR`, `vTheta`, `vPhi` $=dx/dt$), no
  cuatro-velocidades. Para pasar de fórmulas en $\tau$: $v^i=u^i/u^t$.
- Escalas para poner unidades: `rg = bh.r_g`, `T0 = (rg / c).to(u.s)`. Entonces
  `R = r*rg`, `vR = v_r*rg/T0`, `vPhi = Omega*u.rad/T0`.

## 2. Validación 1: ISCO (prógrado y retrógrado)

### Fórmulas (Bardeen, Press & Teukolsky 1972, §II)

Órbitas circulares ecuatoriales de radio $r$:

$$E=\frac{r^{3/2}-2r^{1/2}\pm a}{r^{3/4}\sqrt{r^{3/2}-3r^{1/2}\pm2a}},\qquad
L_z=\pm\frac{r^2\mp2ar^{1/2}+a^2}{r^{3/4}\sqrt{r^{3/2}-3r^{1/2}\pm2a}},\qquad
\Omega\equiv\frac{d\phi}{dt}=\frac{\pm1}{r^{3/2}\pm a}$$

ISCO:

$$Z_1=1+(1-a^2)^{1/3}\left[(1+a)^{1/3}+(1-a)^{1/3}\right],\quad Z_2=\sqrt{3a^2+Z_1^2},\quad
r_{\rm ISCO}=3+Z_2\mp\sqrt{(3-Z_1)(3+Z_1+2Z_2)}$$

Marginalmente ligada: $r_{\rm mb}=2\mp a+2\sqrt{1\mp a}$. Las circulares existen para
$r>r_{\rm ph}$, son ligadas para $r>r_{\rm mb}$ y estables para $r>r_{\rm ISCO}$.

En BPT son las ecs. (2.12) $E$, (2.13) $L$, (2.16) $\Omega$, (2.19) $r_{\rm mb}$, (2.21) $r_{\rm ms}$
**[memoria: confirmar la numeración en el PDF]**.

En el ISCO (Mummery & Balbus 2022, ecs. 7–8) **[fuente]**:

$$E_{\rm ISCO}=\sqrt{1-\frac{2}{3r_I}},\qquad L_{\rm ISCO}=\pm2\sqrt3\left(1\mp\frac{2a}{3\sqrt{r_I}}\right)$$

Frecuencias epicíclicas vistas desde infinito (Okazaki, Kato & Fukue 1987; Kato 1990) **[memoria]**:

$$\Omega_r^2=\Omega^2\left(1-\frac6r\pm\frac{8a}{r^{3/2}}-\frac{3a^2}{r^2}\right),\qquad
\Omega_\theta^2=\Omega^2\left(1\mp\frac{4a}{r^{3/2}}+\frac{3a^2}{r^2}\right)$$

$\Omega_r^2=0$ es exactamente la condición del ISCO.

Caída desde el ISCO, solución cerrada (Mummery & Balbus 2022, ecs. 13 y 18) **[fuente]**. Con
$E=E_{\rm ISCO}$, $L_z=L_{\rm ISCO}$ y $x=r/r_I$:

$$\frac{dr}{d\tau}=-\sqrt{\frac{2}{3r_I}}\left(\frac{r_I}{r}-1\right)^{3/2},\qquad
\tau(r)=\sqrt{\frac{3r_I^3}{2}}\left[3\arcsin\sqrt{x}+(x-3)\sqrt{\frac{x}{1-x}}\right]$$

Solo depende de $a$ a través de $r_I$. Para $a=0$: $r/6=1/(1+12/\phi^2)$ (ec. 16).

### Condiciones iniciales (probadas) **[num]**

```python
bh = rp.Kerr(mass=4.3e6*u.Msun, spin=0.9)
rg = bh.r_g; T0 = (rg / c.c).to(u.s)          # astropy.constants as c
s = +1                                        # +1 prógrada, -1 retrógrada
Omega = s / (r**1.5 + s*bh.spin)              # r en unidades de M
orbit = bh.orbit(R=r*rg, Theta=90*u.deg, Phi=0*u.rad, vPhi=Omega*u.rad/T0)
```

Resultados con `spin=0.9`:

| | prógrada | retrógrada |
|---|---|---|
| $r_{\rm ISCO}$ fórmula | 2.320883041761887 | 8.717352279606489 |
| `bh.r_isco_prograde/retrograde` | 2.3208830417618866 | 8.717352279606489 |
| $E_{\rm ISCO}$ | 0.844247008005536 | 0.961001654354718 |
| $L_{\rm ISCO}$ | 2.099784756123846 | −4.168064196332359 |

- En $r=1.5\,r_{\rm ISCO}$, `get_E()` y `get_Lz()` coinciden con BPT a ~$10^{-15}$.
- Integrando 20 periodos, la deriva relativa de $R$ es ~$10^{-15}$ (`status == 0`).

### Qué mostrar

1. `bh.r_isco_*` contra la fórmula, barriendo $a\in[0,1)$. Límites: $6$ en $a=0$; $1$ y $9$ en $a=1$.
2. Tres órbitas por rama, como la Fig. 7 de PyGRO: estable ($r>r_{\rm ISCO}$), en el ISCO, e
   inestable ($r_{\rm mb}<r<r_{\rm ISCO}$, cae tras acumular error numérico).
3. Cuantitativo, mejor que "se ve circular": perturbar $r\to r(1+10^{-6})$ y medir la frecuencia
   radial contra $\Omega_r$; debe tender a cero al acercarse al ISCO.
4. Analítico de convergencia: arrancar en $r=r_I(1-\epsilon)$ con $E_{\rm ISCO}$, $L_{\rm ISCO}$
   y comparar $\tau(r)$ contra Mummery–Balbus al variar la tolerancia.

Pendiente: para el punto 4 hay que dar $v^r=u^r/u^t$ y $v^\phi=u^\phi/u^t$ con las ecuaciones de §1;
no se probó todavía.

## 3. Validación 2: arrastre de marco

Partícula ecuatorial que cae desde el reposo en infinito: $E=1$, $L_z=0$, $Q=0$. No tiene
momento angular y aun así gira: todo $d\phi/dt$ es arrastre.

### Fórmulas

Salen de sustituir $E=1,\ L_z=0,\ Q=0,\ \theta=\pi/2$ en las ecuaciones de §1 (Carter 1968).
La velocidad angular coincide con la del observador de momento angular cero (ZAMO / LNRF) de
BPT 1972, $\omega=-g_{t\phi}/g_{\phi\phi}$ (ec. 2.? de BPT **[memoria]**):

$$\left(\frac{dr}{d\tau}\right)^2=\frac{2(r^2+a^2)}{r^3},\qquad
\frac{d\phi}{d\tau}=\frac{2a}{r\Delta},\qquad
\frac{dt}{d\tau}=\frac{(r^2+a^2)^2-a^2\Delta}{r^2\Delta}$$

$$\frac{d\phi}{dt}=\omega(r)=\frac{2a}{r^3+a^2r+2a^2}\ \xrightarrow{\ r\to r_+\ }\ \Omega_H=\frac{a}{2r_+}$$

Ángulo acumulado, por cuadratura (diverge logarítmicamente en $r_+$):

$$\phi(r)=\int_r^{r_0}\frac{2a}{r'\Delta(r')}\sqrt{\frac{r'^3}{2(r'^2+a^2)}}\;dr'$$

Límite de campo débil: $\omega\to2a/r^3$ (Lense–Thirring). Para $a=0$: $\phi\equiv0$.

### Condiciones iniciales (probadas) **[num]**

No se puede arrancar en infinito; se arranca en $r_0$ finito ya con la velocidad que tendría:

```python
a, r0 = bh.spin, 1000.0
D  = r0**2 - 2*r0 + a**2
ut = ((r0**2 + a**2)**2 - a**2*D) / (r0**2 * D)
ur = -np.sqrt(2*(r0**2 + a**2) / r0**3)
uph = 2*a / (r0*D)
orbit = bh.orbit(R=r0*rg, Theta=90*u.deg, Phi=0*u.rad,
                 vR=ur/ut*rg/T0, vPhi=uph/ut*u.rad/T0)
sol = orbit.solve(tau_span=(0*T0, 5e4*T0))
```

Resultados con `spin=0.9`, $r_0=1000$:

- Constantes iniciales: $E=1.0$, $L_z=4\times10^{-19}$, $Q=7\times10^{-70}$.
- $|v^\phi/\omega(r)-1|\le2.7\times10^{-9}$ en toda la caída.
- $|u^r/u^r_{\rm teo}-1|\le1.3\times10^{-7}$.
- Último punto: $r=1.4358898963$ contra $r_+=1.4358898944$; $d\phi/dt=0.3133945024$ contra
  $\Omega_H=0.3133945031$.
- $\phi$ total $=20.25240728$ contra cuadratura $20.25240736$ (la cuadratura con `scipy.quad`
  dio aviso de subdivisiones por la singularidad en $r_+$; partir el intervalo o cortar antes).

**Observación para revisar en código:** la integración terminó con `status == -1`
("step size underflow; native status 5") a $r-r_+\approx2\times10^{-9}$, no con `status == 1`
(cruce de horizonte confirmado). Para el notebook conviene cortar en, por ejemplo,
$r=r_+(1+10^{-3})$, o revisar si el evento de horizonte debería dispararse aquí.

### Qué mostrar

1. $d\phi/dt$ contra $r$: puntos integrados sobre la curva $\omega(r)$, con $\Omega_H$ y $2a/r^3$.
2. Trayectoria en el plano $xy$ para varios $a$ (0, 0.5, 0.9, 0.998): la espiral aparece solo con spin.
3. Residuo relativo contra tolerancia (convergencia), como la Fig. 14 de PyGRO.

Extensión opcional, misma idea en una dimensión: caída por el eje ($\theta=0$),
$\dot r^2=E^2-1+2r/(r^2+a^2)$, $\dot t=E(r^2+a^2)/\Delta$. No se probó: relatipy puede rechazar
$\theta=0$ por la singularidad de coordenadas (`get_E` lanza `ValueError` ahí).

## 4. Validación 3: precesión exacta (hecha)

### Teoría

Frecuencias en tiempo de Mino ($d\tau=\Sigma\,d\lambda$; Mino 2003). Con periodos
$\Lambda_r=2\int_{r_p}^{r_a}dr/\sqrt{R}$ y $\Lambda_\theta=4\int_0^{\sqrt{z_-}}d\cos\theta/\sqrt{\Theta}$
(Fujita & Hikida 2009, ecs. 9 y 11 **[fuente]**):

$$\Upsilon_r=\frac{2\pi}{\Lambda_r},\qquad\Upsilon_\theta=\frac{2\pi}{\Lambda_\theta},\qquad
\Upsilon_\phi=\left\langle\frac{a}{\Delta}[E(r^2+a^2)-aL_z]\right\rangle_r-aE+\left\langle\frac{L_z}{\sin^2\theta}\right\rangle_\theta$$

(ec. 13 de Fujita & Hikida; los promedios en $\lambda$ son de Drasco & Hughes 2004). Entonces:

$$\Delta\phi_{\rm peri}=2\pi\left(\frac{\Upsilon_\phi}{\Upsilon_r}-1\right)\ \text{por periodo radial},\qquad
\Delta\phi_{\rm nodo}=2\pi\left(\frac{\Upsilon_\phi}{\Upsilon_\theta}-1\right)\ \text{por periodo polar}$$

Son exactas (sin expansión en $M/p$ ni en $a$), dentro del modelo de partícula de prueba.

Formas cerradas con integrales elípticas (Fujita & Hikida 2009) **[memoria: confirmar en el paper]**,
con raíces $r_1>r_2>r_3>r_4$ de $R$ y $z_-<z_+$ de $\Theta$ en $z=\cos^2\theta$:

$$\Upsilon_r=\frac{\pi\sqrt{(1-E^2)(r_1-r_3)(r_2-r_4)}}{2K(k_r)},\ k_r^2=\frac{(r_1-r_2)(r_3-r_4)}{(r_1-r_3)(r_2-r_4)};\qquad
\Upsilon_\theta=\frac{\pi\sqrt{a^2(1-E^2)z_+}}{2K(k_\theta)},\ k_\theta^2=\frac{z_-}{z_+}$$

Schwarzschild, forma cerrada (parametrización $r=p/(1+e\cos\chi)$ de Cutler, Kennefick & Poisson
1994; también Chandrasekhar 1983, cap. 3) **[num]**:

$$\Delta\phi=4\sqrt{\frac{p}{p-6+2e}}\;K\!\left(\frac{4e}{p-6+2e}\right)-2\pi$$

($K(m)$ con la convención de `scipy.special.ellipk`, parámetro $m=k^2$.)

Campo débil, por órbita (Will 2008, citado como en esta conversación) **[memoria]**:
$A_S=6\pi/p$, $A_J=4\pi\chi/p^{3/2}$, $A_Q=3\pi\chi^2/p^2$;
$\Delta\varpi=A_S-2A_J\cos i-\tfrac12A_Q(1-3\cos^2i)$, $\Delta\Omega=A_J-A_Q\cos i$.
Ojo: $\Delta\varpi$ de Will no es lo mismo que $\Delta\phi_{\rm peri}$ de arriba en órbitas
inclinadas; el avance del pericentro dentro del plano es $2\pi(\Upsilon_\theta/\Upsilon_r-1)$.
Ese mapeo no se hizo; en el notebook la 1PN solo se compara en Schwarzschild.

### Implementación en el notebook

- `exact_precession(spin, E, Lz, Q)`: saca las raíces de $R$ con `np.roots` y evalúa los
  promedios por regla del punto medio con $r=\tfrac{r_a+r_p}2+\tfrac{r_a-r_p}2\cos\xi$ y
  $\cos\theta=\sqrt{z_-}\cos\psi$ (integrandos suaves y periódicos, converge a precisión de
  máquina con `n=4096`). Sin dependencias nuevas.
- `crossing_phases(sol, which)`: encierra cada pericentro (cambio de signo de `uR`) o nodo
  ascendente (de $\pi/2-\Theta$) entre dos muestras e interpola $\Phi$ linealmente.

Comprobaciones de la predicción **[num]**:

- Schwarzschild: cuadratura contra la forma cerrada, coincide a ~$10^{-15}$.
- Kerr ($a=0.9$, $p=12$, $e=0.5$, $x=0.5$): contra `kerrgeopy.StableOrbit(0.9, 12, 0.5, 0.5).mino_frequencies()`,
  coincide a ~$10^{-14}$. Valores: $\Upsilon_r=3.0557284918$, $\Upsilon_\theta=3.8970023142$,
  $\Upsilon_\phi=4.0588920591$; $E=0.9702249839$, $L_z=1.9465955210$, $Q=11.4033404319$.

### Resultados

Schwarzschild, órbita del notebook (`a = 100 r_s`, `e = 0.8`, 5 órbitas, `max_step=T/5000`):

- Exacta: 0.266154212907 rad. 1PN ($6\pi/p$): 0.249441036967 rad (−6.3 %).
- Integrada: diferencia relativa entre $6\times10^{-9}$ y $4\times10^{-7}$.

Kerr genérica (`spin=0.9`, `p=12 r_g`, `e=0.5`, `x=0.5`, 300 periodos keplerianos):

- Pericentro: exacta 2.062703739 rad; media de 140 ciclos 2.061804199 (−4.4×10⁻⁴).
- Nodo: exacta 0.261016850 rad; media de 179 ciclos 0.261048178 (+1.2×10⁻⁴).

### Cosas a tener presentes

- **Los elementos clásicos `a`, `e` de relatipy son osculantes, no puntos de retorno.** Con
  `a = 100 r_s`, `e = 0.8` la geodésica va de $r_p=40$ a $r_a=681.888$ ($p=75.567$, $e=0.8892$).
  Por eso `15*T` kepleriano cubre 5 órbitas y no 15. La 1PN con los elementos nominales da 0.2618.
  Para barrer $(p,e)$ con significado relativista, usar la familia ligada `orbit(p=, e=, x=)`.
- **El residuo en Schwarzschild es de la interpolación, no del integrador.** Interpolación lineal:
  ~$10^{-4}$ con pasos por defecto de Radau, ~$10^{-7}$ con `T/5000`, ~$10^{-10}$ con `T/50000`.
  Con bisección sobre `sol.at` (vectorizada) da ~$3\times10^{-6}$ por defecto y ~$10^{-10}$ con `T/5000`.
  `sol.at` cuesta ~25 ms por llamada con 75 000 muestras sin importar el tamaño de la consulta:
  llamarlo con arreglos, nunca escalar dentro de un `brentq`.
- **En Kerr inclinada el valor ciclo a ciclo oscila** (pericentro entre 1.39 y 2.67 rad). La
  exacta es la media y el error del promedio cae como $1/N$. Para bajar de $10^{-4}$: más ciclos,
  o comparar en órbitas ecuatoriales, donde coincide ciclo a ciclo.
- Una solución con ~$10^6$ muestras quedó truncada (1 000 001 puntos): hay un tope de muestras.

### Falta para cerrarla

- Figura tipo Fig. 6 de PyGRO: precesión contra $p$ para varios $(a,e,x)$, con exacta, 1PN y
  residuo relativo.
- Decidir si el notebook se renombra (ya no es solo Schwarzschild).

## 5. Validación 4: MCMC con S2 y S301

relatipy ya tiene `KerrMcmcModel.get_ra_dec_vr` (tutorial `docs/tutorials/06-mcmc-observables.ipynb`).
El modelo usa proyección de ángulo pequeño y retardo de Rømer en línea recta; no traza fotones
ni incluye Shapiro o lente. Es el mismo nivel de aproximación que PyGRO §3.4.

### S2

Elementos y datos:

| | Gillessen et al. 2017 (Tabla 3) | GRAVITY 2020 |
|---|---|---|
| $a$ | 0.1255″ | 125.058 mas |
| $e$ | 0.8839 | 0.884649 |
| $i$ | 134.18° | 134.567° |
| $\Omega$ | 226.94° | 228.171° |
| $\omega$ | 65.51° | 66.263° |
| $t_P$ | 2002.33 | 2018.379 |
| $P$ | 16.00 yr | 16.0455 yr |
| $M$ | $4.35\times10^6\,M_\odot$ | $4.261\times10^6\,M_\odot$ |
| $R_0$ | 8.33 kpc | 8246.7 pc |

- Columna Gillessen: coincide con el vector de la celda "legacy" del tutorial 06 y con lo que
  usa PyGRO; $M$ y $R_0$ **[fuente: PyGRO §3.4]**, el resto **[memoria]**.
- Columna GRAVITY 2020: **[memoria]**, verificar todo en la Tabla E.1 del paper. Ese paper mide
  la precesión de Schwarzschild de S2, $f_{\rm SP}=1.10\pm0.19$ (~12′ por órbita).
- Datos públicos: las tablas electrónicas de Gillessen et al. 2017 (astrometría NACO y
  velocidades radiales). Son los que usa PyGRO en su Fig. 11. La astrometría de GRAVITY no es pública.
- Convención de ángulos: Gillessen et al. 2009. En GRAVITY 2020 los elementos son osculantes en
  un instante concreto **[memoria: creo que el apocentro de 2010.35; confirmar]**, lo cual
  importa porque relatipy también usa elementos osculantes en el estado inicial.

### S301

Descubierta por la colaboración GRAVITY+ (Nature, 19 de agosto de 2026). Datos del paper **[fuente: PDF de ESO]**:

| Parámetro | Solución 1 | Solución 2 |
|---|---|---|
| $a$ [mas] | 83.0 ± 0.4 | 83.0 ± 0.5 |
| $e$ | 0.9825 ± 0.0013 | 0.9816 ± 0.0006 |
| $i$ [°] | 123.278 ± 0.82 | 121.98 ± 0.26 |
| $\Omega$ [°] | 76.0 ± 2.7 | 259.07 ± 0.40 |
| $\omega$ [°] | 294.7 ± 1.9 | 116.36 ± 0.55 |
| $t_P$ [yr] | 2023.119 ± 0.0095 | 2023.117 ± 0.002 |
| $P$ [yr] | 8.69 ± 0.06 | 8.69 ± 0.08 |
| $\chi^2$ | 29.065 | 29.057 |

- Tabla C2 del paper. Hay dos soluciones porque no hay velocidad radial: la orientación 3D queda
  ambigua. 34 grados de libertad.
- Fijos en el ajuste: $M=4.297\times10^6\,M_\odot$, $R_0=8277$ pc, sin offsets del sistema de referencia.
- Modelo del ajuste: el de S2 de GRAVITY 2020, con Rømer y precesión de Schwarzschild 1PN.
  **No incluye Kerr.**
- 19 épocas astrométricas: 2017.2389; 2021.5759; 2023.2604, 2023.4162, 2023.4298, 2023.5063;
  2024.2325, 2024.2380, 2024.2405, 2024.3281, 2024.4073, 2024.4756, 2024.4837, 2024.5438;
  2025.2070, 2025.3602, 2025.4393, 2025.5159, 2025.5978 (Tabla A1).
- Error astrométrico: ~207 µas por dirección (discretización de la imagen). $m_K=19.3\pm0.3$.
- Pericentro a ~136–142 $r_s$, ~12 UA, ~25 000 km/s (8.3–8.5 % de $c$).
- Precesión de Schwarzschild: 1.9° (sol. 1) y 1.8° (sol. 2) por órbita.
- Lense–Thirring: ~0.11° $\chi\cos\beta$ por revolución en el plano; fuera del plano
  $\propto\chi\sin\beta\sin\lambda$, con $\beta$ el ángulo entre el spin y el momento angular
  orbital y $\lambda$ el ángulo de posición de la proyección del spin en el plano orbital.
  Escalas de tiempo: $T_{\rm SP}\approx1.6\times10^3$ yr, $T_{\rm LT}\approx5.8\times10^4$ yr.
  **Tomar las fórmulas exactas del paper**: la transcripción que obtuve tenía los exponentes dudosos.
- Estudio simulado del paper: ajusta $(\chi,\lambda,\beta)$ con datos sintéticos hasta 2035,
  $\chi=1$ y precisión de 100 µas. Es la referencia natural para comparar el MCMC de relatipy.

### Bloqueos y decisiones para esta validación

1. **No encontré tabuladas las posiciones RA/Dec individuales de S301.** El PDF trae las épocas
   pero, en lo que pude leer, no los valores. Revisar Extended Data / Source Data y la sección
   "Data availability" en la página de Nature (no pude abrirla). Sin eso, S301 solo se puede
   hacer con datos sintéticos generados desde la Tabla C2.
2. Con S301 hay que decidir qué solución usar (o ambas) y qué velocidad radial suponer.
3. Mapear $(\chi,\lambda,\beta)$ del paper al `vec` de relatipy, que va en ejes del observador
   $(\mathrm{Dec},\mathrm{RA},\text{alejándose})$.
4. Propuesta de alcance: (a) S2 con datos de Gillessen 2017, recuperar $M$, $R_0$ y elementos
   con spin fijo en 0, como PyGRO Fig. 11; (b) S301 sintética, recuperar el spin inyectado y
   comparar la incertidumbre con el estudio simulado del paper.

## 6. Referencias

Teoría de geodésicas en Kerr:

- Carter, B. 1968, Phys. Rev. 174, 1559. Separación de las ecuaciones; constante de Carter.
- Bardeen, J. M., Press, W. H. & Teukolsky, S. A. 1972, ApJ 178, 347.
  https://ui.adsabs.harvard.edu/abs/1972ApJ...178..347B (doi:10.1086/151796). Órbitas circulares, ISCO, LNRF.
- Wilkins, D. C. 1972, Phys. Rev. D 5, 814. Órbitas ligadas y esféricas. **[memoria]**
- Chandrasekhar, S. 1983, *The Mathematical Theory of Black Holes* (Oxford). Caps. 3 y 7.
- Schmidt, W. 2002, Class. Quantum Grav. 19, 2743, arXiv:gr-qc/0202090. Frecuencias fundamentales. **[memoria]**
- Mino, Y. 2003, Phys. Rev. D 67, 084027, arXiv:gr-qc/0302075. Tiempo de Mino. **[memoria]**
- Drasco, S. & Hughes, S. A. 2004, Phys. Rev. D 69, 044015, arXiv:astro-ph/0308479. Promedios orbitales. **[memoria]**
- Fujita, R. & Hikida, W. 2009, Class. Quantum Grav. 26, 135002, arXiv:0906.1420.
  https://arxiv.org/abs/0906.1420. Soluciones analíticas de órbitas ligadas.
- Mummery, A. & Balbus, S. 2022, Phys. Rev. Lett. 129, 161101, arXiv:2209.03579.
  https://arxiv.org/abs/2209.03579. Caída desde el ISCO. (Revista **[memoria]**; ecuaciones **[fuente]**.)
- Okazaki, A. T., Kato, S. & Fukue, J. 1987, PASJ 39, 457. Frecuencias epicíclicas en Kerr. **[memoria]**
- Cutler, C., Kennefick, D. & Poisson, E. 1994, Phys. Rev. D 50, 3816. Parametrización $(p,e)$ en Schwarzschild. **[memoria]**

Campo débil y estrellas S:

- Will, C. M. 2008, ApJ 674, L25, arXiv:0711.1677. Precesiones por spin y cuadrupolo. **[memoria]**
- Poisson, E. & Will, C. M. 2014, *Gravity* (Cambridge). Fórmula 1PN (la cita PyGRO).
- Merritt, D., Alexander, T., Mikkola, S. & Will, C. M. 2010, Phys. Rev. D 81, 062002. **[memoria]**
- Gillessen, S. et al. 2009, ApJ 692, 1075. Convención de ángulos orbitales. **[memoria]**
- Gillessen, S. et al. 2017, ApJ 837, 30, arXiv:1611.09144. Datos y elementos de las estrellas S.
- GRAVITY Collaboration 2018, A&A 615, L15. Corrimiento al rojo de S2.
- Do, T. et al. 2019, Science 365, 664. Corrimiento al rojo de S2 (Keck).
- GRAVITY Collaboration 2020, A&A 636, L5. Precesión de Schwarzschild de S2.
- GRAVITY Collaboration 2022, A&A 657, L12. Distribución de masa; S2, S29, S38, S55. **[memoria]**
- Grould, M. et al. 2017, A&A 608, A60. Efectos relativistas en los observables de S2. **[memoria]**
- GRAVITY+ Collaboration (Abd El Dayem, K. et al.) 2026, Nature, "Discovery of a star sensitive
  to the spin of Sagittarius A*", doi:10.1038/s41586-026-10894-w, bibcode 2026Natur.657..359A.
  PDF: https://www.eso.org/public/archives/releases/sciencepapers/eso2612/eso2612a.pdf
- "Probing departures from Schwarzschild geometry with the S301 star and kilometre-baseline
  interferometry", arXiv:2610.07123. **No lo pude leer** (el servidor rechazó la consulta);
  apareció en la búsqueda y parece directamente relevante.
- Jin et al. 2022, MNRAS 509, 4608. Recuperación bayesiana de órbitas (referencia de método ya usada).

Software:

- Della Monica, R. 2025, PyGRO, arXiv:2504.20152. https://arxiv.org/abs/2504.20152. Modelo de la sección de validación.
- Park, S. & Nasipak, Z. 2024, KerrGeoPy, JOSS, doi:10.21105/joss.06587, arXiv:2406.01413.
  https://joss.theoj.org/papers/10.21105/joss.06587
- Foreman-Mackey, D. et al. 2013, PASP 125, 306. emcee.

## 7. Pendientes de verificación

- Numeración de ecuaciones en BPT 1972 y fórmulas elípticas de Fujita & Hikida (§2 y §4).
- Toda la columna GRAVITY 2020 de S2 y el instante de osculación.
- Fórmulas de Lense–Thirring del paper de S301, con exponentes y $c$ bien puestos.
- Disponibilidad de la astrometría de S301.
- Todas las entradas marcadas **[memoria]** en §6: volumen, página y año.
- Compatibilidad de los notebooks con la versión de relatipy que esté en PyPI al publicar.
- **Dependencias de PyGRO para los pares (por revisar).** `pygro_engine` usa
  `backend="autowrap"` (el de la doc de PyGRO): sympy genera C y lo compila con Cython, así que
  el revisor necesita compilador C + `cython` + `setuptools` (falló en `thesis/.venv`, Python
  3.14, por falta de Cython). Alternativa: `backend="lambdify"`, Python puro, sin compilar,
  más lento (no medido). Decidir: dejar autowrap y documentar requisitos, cambiar a lambdify,
  o probar autowrap y caer a lambdify si falla (PyGRO solo cae solo con `OSError`, no con
  `CodeWrapError`). Si cambia el backend, re-medir tiempos; la precisión no debería cambiar
  (mismas ecuaciones), pero verificar.
