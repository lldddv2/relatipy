# Boyer-Lindquist to Cartesian coordinates

$$
\begin{align*}
\rho &= \sqrt{r^2 + a^2} \\
x &= \rho \sin\theta \cos\phi \\
y &= \rho \sin\theta \sin\phi \\
z &= r \cos\theta
\end{align*}
$$

en nuestra notación:

$$
\begin{align*}
\rho &= \sqrt{(x^1)^2 + a^2} \\
x'^1 &= \rho \sin x^2 \cos x^3 \\
x'^2 &= \rho \sin x^2 \sin x^3 \\
x'^3 &= x^1 \cos x^2
\end{align*}
$$

# Inverse transformation
Es claro que:
$$
x^2 = \arccos\left(\frac{x'^3}{x^1}\right)
$$
y 
$$
x^3 = \arctan\left(\frac{x'^2}{x'^1}\right)
$$

Para $r$:

$$
R^2 = (x'^1)^2 + (x'^2)^2 + (x'^3)^2
$$

así:
$$
\begin{align*}
R^2 &=  (x'^1)^2 + (x'^2)^2 + (x'^3)^2 \\
&= (\rho \sin x^2 \cos x^3)^2 + (\rho \sin x^2 \sin x^3)^2 + (x^1 \cos x^2)^2 \\
&= \rho^2 \sin^2 x^2 (\cos^2 x^3 + \sin^2 x^3) + (x^1)^2 \cos^2 x^2 \\
&= \rho^2 \sin^2 x^2 + (x^1)^2 \cos^2 x^2 \\
&= ( (x^1)^2 + a^2 ) \sin^2 x^2 + (x^1)^2 \cos^2 x^2 \\
&= (x^1)^2 (\sin^2 x^2 + \cos^2 x^2) + a^2 \sin^2 x^2 \\
&= (x^1)^2 + a^2 \sin^2 x^2
\end{align*}
$$

De modo que:
$$
(x^1)^2 = R^2 - a^2 \sin^2 x^2
$$

Pero como $\sin^2 x^2 = 1 - \cos^2 x^2$ y además $\cos x^2 = \frac{x'^3}{x^1}$, tenemos que:

$$
\begin{align*}
(x^1)^2 &= R^2 - a^2 (1 - \cos^2 x^2) \\
&= R^2 - a^2 + a^2 \frac{(x'^3)^2}{(x^1)^2} \\
&= \frac{(x^1)^2 (R^2 - a^2) + a^2 (x'^3)^2}{(x^1)^2}
\end{align*}
$$

Por lo que:
$$
(x^1)^4 - (R^2 - a^2) (x^1)^2 - a^2 (x'^3)^2 = 0
$$

Las raices para $(x^1)^2$ son:

$$
(x^1)^2 = \frac{(R^2 - a^2) \pm \sqrt{(R^2 - a^2)^2 + 4 a^2 (x'^3)^2}}{2}
$$

Nos interesan raices reales, por lo que debemos garantizar que $(x^1)^2 \geq 0$.

Note que, para cualquier $x'^3 \neq 0$, y $a, R \in \mathbb{R}$, se cumple que:

$$
\begin{align*}
\sqrt{(R^2 - a^2)^2 + 4 a^2 (x'^3)^2} &\geq \sqrt{(R^2 - a^2)^2} \\
&\geq |R^2 - a^2| \\
\Rightarrow 0 &\geq |R^2 - a^2| - \sqrt{(R^2 - a^2)^2 + 4 a^2 (x'^3)^2} \\
\end{align*}
$$


Esto es, no tomamos la raíz negativa. Si $x'^3 = 0$. Para la raiz positiva, se cumple que:

$$
\begin{align*}
\sqrt{(R^2 - a^2)^2 + 4 a^2 (x'^3)^2} &\geq 0 \\
&\geq - \sqrt{(R^2 - a^2)^2} \\
&\geq - |R^2 - a^2| \\
\Rightarrow |R^2 - a^2| + \sqrt{(R^2 - a^2)^2 + 4 a^2 (x'^3)^2} &\geq 0 \\
\end{align*}
$$

Así, siempre que $R^2 \geq a^2$, se cumple que $(x^1)^2 \geq 0$ para la raíz positiva. En geneal, asumiremos que $R^2 \geq a^2$.

Finalmente, la transformación inversa es:
$$
\begin{align*}
R^2 &= (x'^1)^2 + (x'^2)^2 + (x'^3)^2 \\
x^1 &= \sqrt{\frac{(R^2 - a^2) + \sqrt{(R^2 - a^2)^2 + 4 a^2 (x'^3)^2}}{2}} \\
x^2 &= \arccos\left(\frac{x'^3}{x^1}\right) \\
x^3 &= \arctan\left(\frac{x'^2}{x'^1}\right)
\end{align*}
$$