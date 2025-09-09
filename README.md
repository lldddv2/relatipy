# KerrPy 🌀

Herramienta simbólica, numérica y visual para explorar la geometría y dinámica alrededor de **agujeros negros de Kerr**.  

---

## 📑 Tabla de contenidos
- [KerrPy 🌀](#kerrpy-)
  - [📑 Tabla de contenidos](#-tabla-de-contenidos)
  - [🌌 Introducción](#-introducción)
  - [🏗️ Arquitectura del proyecto](#️-arquitectura-del-proyecto)
  - [Inicializar proyecto](#inicializar-proyecto)

---

## 🌌 Introducción
KerrPy tiene como objetivo brindar un entorno unificado para estudiar **agujeros negros de Kerr** combinando:
- **Cálculos simbólicos** (tensores, invariantes, geodésicas).
- **Simulaciones numéricas** (órbitas, trayectorias de partículas y fotones).
- **Visualización interactiva** (potenciales efectivos, mapas de curvatura, trayectorias).  

Este proyecto busca servir como herramienta de **investigación**.

---

## 🏗️ Arquitectura del proyecto
El proyecto se organiza en módulos separados:
```
/docs          → documentación teórica y notas
/symbolic      → cálculos simbólicos 
/numeric       → cálculos numéricos 
/visualization → visualización y análisis interactivo 
```

## Inicializar proyecto
Garantizando que tiene `docker-compose` instalado, puede iniciar el proyecto con:

```bash
docker-compose build
```

Para la parte simbólica:
```
docker-compose run symbolic
```

Para la parte numérica:
```
docker-compose run numeric
```

Para la parte de visualización:
```
docker-compose up visualization
```