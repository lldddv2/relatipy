# 1. Convención de ramas.
La rama principal del proyecto es `main`. La rama secundaria es `dev`, toda nueva funcionalidad o corrección de errores debe hacerse en ramas derivadas de `dev`. 

Se dividen en tres tipos de ramas:
- **feature**: para nuevas funcionalidades, deben derivar de `dev` y al terminarse se hace merge a `dev`.
- **refactor**: para refactorizaciones de código, deben derivar de `dev` y al terminarse se hace merge a `dev`.
- **bugfix**: para corrección de errores, deben derivar de `dev`
- **docupdate**: para actualizaciones de documentación, deben derivar de `dev` y al terminarse se hace merge a `dev`.

Cada rama tiene un código de la siguiente forma:`type/Classification-ShortDescription-###`, los `type`s on los tipos de ramas mencionados antes, `Classification` sólo puede ser:
- `NUM`: para la parte numérica del proyecto.
- `SYM`: para la parte simbólica del proyecto.
- `VIS`: para la parte de visualización del proyecto.
- `DOC`: para la parte de documentación del proyecto.
`ShortDescription` es una breve descripción de la funcionalidad esto es, no más de 5 palabras, y `###` representa el número asociado a la cantidad de funcionalidad asociada a un mismo tema.

Por ejemplo, si creamos una nueva funcionalidad numérica de los símbolos de Christoffel, la rama se llamaría: `feature/NUM-ChristoffelSymbols-001`.

# 2. Convención de commits.
Los commits deben empezar con alguna de las siguientes etiquetas:
- `[FEAT]`: para nuevas funcionalidades.
- `[FIX]`: para corrección de errores.
- `[DOC]`: para actualizaciones de documentación.
- `[REFACTOR]`: para refactorizaciones de código.
- `[TEST]`: para adición o modificación de pruebas.
- `[STYLE]`: para cambios de estilo o formato (ej: formateo con `black`).

Además, los mensajes de commit deben ser claros y concisos, describiendo brevemente el cambio realizado. No se debe entrar en detalles teóricos. 

- Los commits deben ser atómicos, es decir, cada commit debe representar un único cambio o una única funcionalidad.