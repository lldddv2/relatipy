# Traspaso: S301 y amigos con Kerr exacto y potencial no esférico

Estado al 2026-10-25. Este documento basta para continuar tras compactar la
conversación.

## 0. Preferencias del usuario (obligatorias)
- Responder en español. Notas internas en estilo "caveman wenyan" y respuestas
  finales en "caveman ultra", conservando números y citas.
- Nunca usar Sonnet. Haiku (esfuerzo máximo) para búsquedas y extracción; Opus
  permitido para juicio y cálculo. Delegar en subagentes: el contexto principal
  es caro.
- No modificar el código de la librería sin permiso. Trabajar en `research/`.
- Hacer commit y push al branch `claude/funny-galileo-7lm47r` (el hook de stop
  lo exige). No abrir PR.
- arxiv.org está bloqueado en el entorno: leer papers con el conector alphaXiv
  (`mcp__alphaXiv__answer_pdf_queries`, cargar con ToolSearch). Las búsquedas
  van con WebSearch.
- El usuario quiere discutir a fondo el plan antes de ejecutar más.

## 1. Dirección elegida
**Ajuste conjunto de órbita completa (sin promediar) con Kerr exacto más un
potencial newtoniano no esférico (disco delgado o cúspide achatada) para S301 y
sus "amigos" S2, S55, S38 y S29, por MCMC.** El objetivo es medir el vector de
espín de Sgr A* separándolo del fondo de masa extendida.

Por qué hace falta la librería:
- Kerr más un potencial no es integrable, así que las soluciones analíticas
  (Fujita–Hikida, KerrGeoPy) no sirven y hay que integrar numéricamente.
- El MCMC tiene unos 45 parámetros con 5 estrellas, así que la velocidad es
  necesaria.

## 2. Base: Piran et al. 2026, "S301 and friends" (arXiv:2607.24931), leído completo
- Autores: Piran, Amaro-Seoane, Aytac, Bourdarot, Burkert, Calderón, Cuadra,
  Eisenhauer, Genzel, Gillessen, Joharle, Mang, Naab, Ott, Perets, Ribeiro,
  Sadun Bordoni, Sari, Thiel. Es decir, el núcleo de GRAVITY.
- Método: torque newtoniano **secular y promediado por órbita** (analítico, más
  validación con anillos de disco promediados con peso kepleriano y
  ablandamiento H = 0.01 r'). **No** integran órbitas completas, **no** ajustan
  datos simulados, **no** hacen Fisher ni MCMC y **no** dan precisión de espín.
- Física clave: en órbitas excéntricas el torque nodal del disco lo domina el
  **apocentro** y el Lense–Thirring el **pericentro** (∝ (r_g/r_p)^{3/2}).
  - S301, S2, S55 y S38 tienen apocentros parecidos (31.7–50.4 mil r_g), así que
    sienten un torque de disco comparable.
  - El LT de S2 es ~30× menor que el de S301.
  - Por tanto S2, S55 y S38 calibran el disco y S29 fija el perfil radial. Se
    necesitan las cuatro, porque sus planos desalineados triangulan la
    orientación del disco.
- Disco fiducial: razor-thin, Σ ∝ R^α con α = −1, borde exterior en el apocentro
  de S301, M_disk = 10³ Msun (cota superior; el disco maximiza el cuadrupolo),
  χ = 1.
- Tabla 1 (S301 / S2), cociente |ΔΩ_disco/ΔΩ_LT|:

  | i | ω | S301 | S2 | ΔΩ_disco(S2)/ΔΩ_disco(S301) |
  |---|---|---|---|---|
  | 10° | 0° | 0.41 | 13 | 1.0 |
  | 10° | 90° | 2.5 | 19 | 0.26 |
  | 30° | 0° | 0.36 | 8.9 | 0.85 |
  | 30° | 90° | 1.4 | 14 | 0.34 |
  | 60° | 0° | 0.19 | 4.2 | 0.75 |
  | 60° | 90° | 0.61 | 6.6 | 0.37 |

- Valores para χ = 1 y orientación favorable: ΔΩ_LT(S301) = 9.6e-4 rad/órbita;
  ΔΩ_LT(S2) = 3.3e-5.
- Presupuesto nodal de S301: LT 9.6e-4; disco ~2.8e-3, con rango de 2e-4 a
  1.4e-2 según la geometría.
- Presupuesto apsidal de S301: Schwarzschild 3.4e-2; LT 2.0e-3; disco 7e-4
  (variable en el tiempo); 2PN 2.9e-4, que debe restarse a precisión 2PN;
  cuadrupolo de espín ~4 % del LT, despreciable.
- El plano orbital (nodal + deriva de inclinación) mide χ sinβ y λ "en pocas
  órbitas". La componente χ cosβ sale del apsidal y tarda décadas.
- El discriminante temporal (la rotación apsidal de S301, 1.95°/órbita, cambia
  el torque del disco) actúa sobre "a few dozen orbital periods".
- La no detección de nodal en S2 implica ≲1.7e-4 rad/órbita, lo que reduce ~3×
  el cuadrupolo máximo supuesto.
- Huecos que **ellos mismos** declaran:
  1. "this residual geometric mismatch must be fitted jointly within the global
     orbital model" (no lo hacen).
  2. Granularidad en el pericentro de S301: "we leave this to future work". Se
     apoyan en Sadun Bordoni+2025: con perturbadores de 100 Msun la masa
     encerrada puede sesgarse hasta ×6, y el retroceso browniano de Sgr A* llega
     a 6 μas.
  3. Ambigüedad entre las dos soluciones orbitales de S301 (A/B), que difieren
     en orientación.
- Tabla 3 (a [r_g], e, r_a): S301 16010, 0.9825, 31740; S55 20334, 0.7267, 35110;
  S2 24100, 0.8844, 45400; S38 27763, 0.8145, 50386; S29 77421, 0.9693, 152465.
  Con M = 4.297e6 Msun y R0 = 8.28 kpc, 1″ ≈ 195 350 r_g.

## 3. Qué aporta el proyecto frente a Piran+
1. El **ajuste conjunto global** que piden y no hacen: órbita completa, MCMC con
   las 5 estrellas y el disco marginalizado (masa, normal y perfil).
2. **Información dentro de cada órbita**, que el promedio secular descarta: el LT
   gira el plano en el pericentro (casi impulsivo) y el disco gradualmente en la
   mitad del apocentro. Podría separarlos en 1–2 órbitas, no en décadas.
   Antecedente: Heißel+2022 hizo lo análogo en el plano para S2.
3. **RV** (ELT/MICADO ~1 km/s). En mi Fisher anterior de S301 la mayor parte de
   la información del espín venía de la RV.
4. Cualquier geometría: disco horario, cúspide achatada, ley de potencias libre.
5. Más adelante, perturbadores discretos (granularidad) dentro de la
   integración.

Riesgo: GRAVITY tiene los datos y hará el ajuste global (probablemente con su
código PN). Ventajas propias: código abierto, Kerr exacto y ser el primer
forecast cuantitativo. Calendario: pericentros de S301 en 2031 y ~2040.

## 4. Otros papers de la misma ola
- **Descubrimiento de S301**: GRAVITY 2026, Nature, DOI
  10.1038/s41586-026-10894-w, arXiv:2607.12664. K ≈ 19.3, P ≈ 8.7 años,
  e ≈ 0.982, r_p ≈ 280 r_g, ~25 600 km/s en el pericentro; soluciones A/B.
- **Disco en el Centro Galáctico** (arXiv:2607.26713, Foschi, Vincent, Paumard,
  Perrin, según una fuente secundaria): aceleración analítica de un disco
  delgado con Σ en ley de potencias; elementos osculantes; marco newtoniano más
  precesión de Schwarzschild PN; disco alineado con el disco estelar horario.
  Para S301 con χ = 1 y θ = π/2, el disco es ~9× menor que el LT, pero la masa es
  incierta.
- **Campaña de dos pasos** (arXiv:2609.06936): el frame-dragging contamina la
  precesión de Schwarzschild en ~5 % (χ = 0.90 ± 0.06, del método de outflow).
- **Rotación estelar de S301** como giroscopio (arXiv:2607.26134, 2609.06977).
- **GRAVITY 2024** (A&A 692 A242, arXiv:2409.12261): masa extendida esférica
  < 1200 Msun (1σ) dentro de S2; f_SP(S2) ≈ 0.91–0.92 ± 0.13.
- **Merritt+2010** (PRD 81 062002): una cúspide achatada con
  M_enc/M ~ (r_g/a)^{3/2} da nodal comparable al LT. Iorio 2011 calculó discos.
  Will 2008 propuso usar estrellas con distinto r_p.

## 5. Números propios (en `results/`)
Precesión por órbita en rad. Masa de la cúspide normalizada al límite de GRAVITY
o al cúmulo nuclear:

| | S2 | S301 |
|---|---|---|
| 1PN | 3.5e-3 | 3.4e-2 |
| 2PN | 3.1e-6 | 2.7e-4 |
| LT apsidal (a*=0.7) | 6.8e-5 | 2.0e-3 |
| LT nodal (a*=0.7) | 2.3e-5 | 6.6e-4 |
| cúspide esférica (retrógrada, en el plano) | −0.6…−6.9e-4 | −0.1…−1.9e-4 |

- La Galaxia a gran escala es despreciable: marea de 1e-12 a 1e-8.
- Fisher de S301 (a* = 0.7, 60 ejes, ángulos sustitutos, **sin disco**):
  - 2026–2035, 100 μas / 1 km/s: σ(a*) mediana 0.60.
  - 2026–2040, 100 μas / 1 km/s: 0.19.
  - 2026–2040, 10 μas / 1 km/s: 0.079.
  - La RV domina la información.
- Lección de MCMC: con 40 walkers en 17D el ensemble colapsó. Usar ≥4 walkers
  por dimensión e inicializar desde Fisher.

## 6. Resultado lateral, pendiente de verificar
Truth = Kerr exacto + lente + Shapiro a primer orden; S2 con datos actuales sin
ruido (`round2/taskA_*`). Un modelo tipo GRAVITY (1PN armónico × f_SP, recta,
sin luz) sesga f_SP en **−0.15**, unas 1.25 veces la σ publicada; ocurre igual
con a* = 0 y con 0.7. El desglose atribuye casi todo a la luz omitida: sin lente
−0.096; sin Shapiro ni deflexión en la RV −0.081; sin ninguno de los tres
−0.027. Va en la dirección del 0.91 medido por GRAVITY. **Falta comprobar si el
modelo real de GRAVITY ya incluye lente y Shapiro** (GRAVITY 2024 / 2020).

## 7. Librería: lo relevante para este proyecto
- Falta una **fuerza externa**: añadir −∇Φ proyectado ortogonal a u^μ en el RHS
  de C. Con Φ_ext/Φ_BH ~ 1e-4 y r ≫ r_g el acoplamiento es inequívoco.
- Con potencial no esférico **no se conservan Q, y E/Lz tampoco si no es
  axisimétrico respecto al espín**. `projection_radau` habría que adaptarlo o
  no usarlo.
- `dop853` es el mejor integrador por defecto (`radau` es ~12× más lento y
  disipativo). `projection_radau` falla en el punto de retorno radial con
  rtol ≥ 1e-10.
- Observables: proyección recta en BL y redshift `1+z = u^t + u_away` sin
  factores métricos (`native/src/geodesic/physic/kerr_observables.c:143`). Error
  de RV (M/r)·v_los·tan²(ψ/2), entre 3 y 609 km/s en el pericentro de S301.
  Para datos reales a 1 km/s hace falta una función T de primer orden (Rømer +
  Shapiro; ≤1.5 km/s con ψ ≤ 120°) o rayos exactos. El núcleo C ya tiene
  `rp_kerr_null_geodesic_rhs`.

## 8. Trabajo en curso (subagentes) y archivos
- **Prototipo Kerr+disco (Opus)** escribe en `research/purposes_2026-10/disk/`.
  Dinámica PN validada (1PN + espín-órbita, 2PN omitido igual en truth y modelo)
  con aceleración del disco por suma de anillos tabulada. Preguntas:
  - Q1: σ del vector de espín (cartesiano; χ⊥ y χ∥) para S301 sola sin disco,
    S301 sola con disco marginalizado, y las 5 estrellas con disco
    marginalizado; con y sin RV (1 y 5 km/s); con 1 o 2 pericentros.
  - Q2: información dentro de la órbita frente a secular.
  - Q3: sesgo de espín si se ignora el disco.
- **Extracción de papers (haiku)** escribe en
  `/tmp/claude-0/-home-user-relatipy/9021de26-e198-5a73-ad04-6e893f8ac1b1/scratchpad/extract_papers.md`:
  elementos de S301 A/B, fórmula del disco de Foschi+, elementos de S2, S29, S38
  y S55 (GRAVITY 2024), precisiones y campaña.
- Scripts validados en el scratchpad `work/`, con copias en
  `research/purposes_2026-10/scripts/`: `pn_fast.py` (numba PN + Rømer + Hermite),
  `pn.py`, `common.py`, `fisher.py`, `chk.py` (rayos exactos de Schwarzschild).
- Informe general previo: `research/purposes_2026-10/REPORT.md` (otros
  propósitos ya descartados).
