# 1. Convenciones generales
Este documento establece las convenciones y notaciones utilizadas en el proyecto **KerrPy** para asegurar consistencia y legibilidad en el código y la documentación.

El idioma principal del proyecto es el inglés, se permite el uso del español únicamente en documentación teórica.

## 1.1. Nombrar variables y funciones
- **Variables**: Utilizar `snake_case` para nombrar variables (ej: `event_horizon_radius`).
- **Funciones**: Utilizar `snake_case` para nombrar funciones (ej: `calculate_geodesics()`).
- **Clases**: Utilizar `PascalCase` para nombrar clases (ej: `KerrMetric`).
- **Constantes**: Utilizar `UPPER_SNAKE_CASE` para nombrar constantes (ej: `SPEED_OF_LIGHT`).
- **Módulos y paquetes**: Utilizar `snake_case` para nombrar módulos y paquetes (ej: `symbolic_calculations`).

**Nota:** Se debe establecer el tipo de variable (int, float, str, list, dict, etc.) en la declaración de la variable o función.

### 1.1.1. Casos particulares
- **Variables indiciales**: Utilizaremos letras griegas y latinas según la convención de relatividad general **ÚNICAMENTE** para variables que representen índices, está prohibido su uso para cualquier otra variable como iteradores, contadores, fuera de este contexto.
- **Variables reservadas**: En python existen ciertas variables reservadas que pueden entrar en conflicto con nombres de variables o funciones, por ejemplo `lambda`. Siempre que sea posible, se debe evitar su uso, en caso de se necesario su uso explicito, pueden usarse variaciones como `lambda_` o `lambda_val`.

## 1.2. Estructura del proyecto
El proyecto se organiza en módulos separados para facilitar el mantenimiento y la escalabilidad:
```
/docs          → documentación teórica y notas
/symbolic      → cálculos simbólicos
/numeric       → cálculos numéricos
/visualization → visualización y análisis interactivo
```

Cada módulo debe contener un archivo `README.md` que describa su propósito y uso, una carpeta `tests/` para pruebas unitarias y de integración, y un archivo `main.py` como punto de entrada.

No se permite la importación cruzada entre módulos, cada módulo debe ser independiente.

## 1.3. Documentación
### 1.3.1. Docstrings
- Todas las funciones, clases y módulos deben incluir docstrings siguiendo el formato de [Google Python Style Guide](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings).
- Todas las docstrings deben estar en inglés para asegurar la accesibilidad internacional del proyecto.
- El uso de ecuaciones debe hacerse mediante LaTeX dentro de las docstrings, utilizando delimitadores adecuados para ecuaciones, ```:math:`expr` ``` para inline equations, o ```.. math::   expr``` para block equations.
- **TODO** docstring debe empezar con `r""" """`, para evitar problemas con caracteres especiales y secuencias de escape.

Aquí ejemplo de docstring para una función:
```python
from typing import Sequence

def compute_friedmann_equation(
    scale_factor: float,
    density: float,
    pressure: float,
    speed_of_light: float,
    gravitational_constant: float,
    cosmological_constant: float,
) -> float:
    r"""Computes the Friedmann acceleration equation.

    This function calculates the second derivative of the scale factor
    :math:`\ddot{a}` normalized by the scale factor :math:`a`, which
    describes the acceleration (or deceleration) of the universe’s
    expansion. The equation is derived from Einstein’s field equations
    under the assumption of a homogeneous and isotropic universe.

    .. math::

        \frac{\ddot{a}}{a} =
        -\frac{4 \pi G}{3}\left(\rho + \frac{3p}{c^2}\right) +
        \frac{\Lambda c^2}{3}

    Args:
      scale_factor:
        The cosmological scale factor :math:`a` (dimensionless).
      density:
        The matter-energy density :math:`\rho` (kg·m⁻³).
      pressure:
        The pressure :math:`p` (Pa).
      speed_of_light:
        The speed of light :math:`c` (m·s⁻¹).
      gravitational_constant:
        Newton’s gravitational constant :math:`G` (m³·kg⁻¹·s⁻²).
      cosmological_constant:
        The cosmological constant :math:`\Lambda` (m⁻²).

    Returns:
      The normalized acceleration :math:`\ddot{a} / a` (s⁻²).

    Raises:
      ValueError: If the scale factor is zero or negative.

    Examples:
      >>> compute_friedmann_equation(
      ...     scale_factor=1.0,
      ...     density=1e-26,
      ...     pressure=0.0,
      ...     speed_of_light=3e8,
      ...     gravitational_constant=6.67430e-11,
      ...     cosmological_constant=1e-52
      ... )
      -2.0963691786580865e-36
    """
    pass
```

### 1.3.2. Comentarios
- Los comentarios deben ser claros y concisos, explicando el "por qué" detrás de decisiones complejas en el código. 
- **NO** deben explicar la física o matemáticas detrás del código, para eso está la documentación teórica en `/docs/relativity`.

## 1.4. Testing
- Todas las funciones y clases deben tener pruebas unitarias que precedan su implementación.
- Utilizar `pytest` como framework de pruebas.

## 1.5. Formateo y linting
- Utilizar `black` para formateo automático del código.
- Utilizar `flake8` para linting y asegurar la calidad del código.
- Utilizar `isort` para ordenar automáticamente las importaciones.


# 2. Preferencias no obigatorias pero recomendadas
- Funciones atómicas: Cada función debe realizar una única tarea o cálculo específico.