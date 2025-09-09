# Imagen base ligera
FROM python:3.11-slim

# Metadata (
LABEL maintainer="Luis Daniel Díaz <luis.diaz@udea.edu.co>" \
      description="Contenedor desarrollar kerrpy" \
      version="1.0"

# Configuración de entorno
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/home/appuser/.local/bin:$PATH"

# Instalar dependencias del sistema
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        git \
    && rm -rf /var/lib/apt/lists/*

# Crear usuario no root
RUN useradd -m appuser

# Directorio de trabajo
WORKDIR /app

# Copiar requirements primero (mejora cacheo de capas)
COPY requirements.txt .

# Instalar dependencias Python como usuario normal
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# Copiar el resto del código
COPY . .

# Cambiar al usuario sin privilegios
USER appuser

# Comando por defecto (sobrescribible en docker-compose)
CMD ["python", "main.py"]
