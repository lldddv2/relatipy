# Kerr black hole metric
EL elemento de linea de Kerr es:

$$
\begin{align*} 
  d s^2=
  & +\left(1-\frac{R_sr}{\Sigma}\right) c^2 d t^2 \\
  &-\frac{\Sigma}{\Delta} d r^2 \\
  &-\Sigma d \theta^2 \\
  &-\left(r^2+a^2+\frac{R_s r a^2}{\Sigma} \sin ^2 \theta\right) \sin ^2 \theta d \phi^2 \\
  &+\frac{2 R_s r a \sin ^2 \theta}{\Sigma} c d t d \phi
\end{align*}
$$

donde $R_s=2GM/c^2$ es el radio de Schwarzschild, $a=J/(Mc)$ es llamado el *parámetro de Kerr* o parámetro de momentum angular ($J$ es el momentum angular total del cuerpo), y las cantidades auxiliares $\Sigma$ y $\Delta$ se pueden escribir como:

$$
\begin{aligned}
    & \Sigma=r^2\left(1+\frac{a^2}{r^2} \cos ^2 \theta\right) \\
    & \Delta=r^2\left(1-\frac{sr+a^2}{r^2}\right)
\end{aligned}
$$

Para hacerlo compatible con nuestra convención, reemplazaremos: 
- $x^0 = ct$
- $x^1 = r$
- $x^2 = \theta$
- $x^3 = \phi$

De este modo:

$$
\begin{align*} 
  d s^2=
  & +\left(1-\frac{R_sx^1}{\Sigma}\right) (d x^0)^2 \\
  &-\frac{\Sigma}{\Delta} (d x^1)^2 \\
  &-\Sigma (d x^2)^2 \\
  &-\left((x^1)^2+a^2+\frac{R_s x^1 a^2}{\Sigma} \sin ^2 (x^2)\right) \sin ^2 (x^2) (d x^3)^2 \\
  &+\frac{2 R_s x^1 a \sin ^2 (x^2)}{\Sigma} d x^0 d x^3
\end{align*}
$$

con:

$$
\begin{align*}
    & \Sigma=(x^1)^2\left(1+\frac{a^2}{(x^1)^2} \cos ^2 (x^2)\right) \\
    & \Delta=(x^1)^2\left(1-\frac{R_sx^1+a^2}{(x^1)^2}\right)
\end{align*}
$$

Para hacerlo más legible diremos:
- $A = 1 - \frac{R_s x^1}{\Sigma}$
- $B = -\frac{\Sigma}{\Delta}$
- $C = -\Sigma$
- $D = -\left((x^1)^2 + a^2 + \frac{R_s x^1 a^2}{\Sigma} \sin^2(x^2)\right) \sin^2(x^2)$
- $E = \frac{2 R_s x^1 a \sin^2(x^2)}{\Sigma}$ 

De modo que la componente de línea queda:
$$
\begin{equation}
d s^2 = A (d x^0)^2 + B (d x^1)^2 + C (d x^2)^2 + D (d x^3)^2 + E d x^0 d x^3
\end{equation}
$$

Y la métrica es:
$$
\begin{equation}
    
g_{\mu \nu} = \begin{pmatrix}
A & 0 & 0 & \frac{E}{2} \\
0 & B & 0 & 0 \\
0 & 0 & C & 0 \\
\frac{E}{2} & 0 & 0 & D
\end{pmatrix}
\end{equation}
$$

Para hallar el tensor de derivadas, podemos usar `sympy` y tener la derivada analítica. La ventaja de esto, es que, por muy larga que sea la expresión, tenerla de forma analítica nos permite portabilidad a otros lenguajes de programación, como `C` o `Fortran`.
