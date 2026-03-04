import plotly.graph_objects as go
import numpy as np

def construct_sphere_surface(R_s, color='black', opacity=1.0):
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x_s = R_s * np.outer(np.cos(u), np.sin(v))
    y_s = R_s * np.outer(np.sin(u), np.sin(v))
    z_s = R_s * np.outer(np.ones(np.size(u)), np.cos(v))

    surface = go.Surface(
        x=x_s,
        y=y_s,
        z=z_s,
        opacity=opacity,  # opacidad completa para esfera negra
        showscale=False,
        name="Schwarzschild radius",
        colorscale=[[0, color], [1, color]]
    )

    return surface

def construct_circular_plane_surface(semi_side, color='black', opacity=1.0):
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x_s = semi_side * np.outer(np.cos(u), np.sin(v))
    y_s = semi_side * np.outer(np.sin(u), np.sin(v))
    z_s = np.zeros_like(x_s)

    surface = go.Surface(
        x=x_s,
        y=y_s,
        z=z_s,
        opacity=opacity,
        colorscale=[[0, color], [1, color]],
        showscale=False,
    )

    return surface

def construct_path(path, color='red', opacity=1.0, label="line"):
    return go.Scatter3d(
        x=path[1],
        y=path[2],
        z=path[3],
        mode='lines',
        line=dict(color=color, width=2),
        opacity=opacity,
        name=label,
    )

def construct_point(x, y, z, color='red', size=10, opacity=1.0, label="point"):
    return go.Scatter3d(
        x=[x],
        y=[y],
        z=[z],
        mode='markers',
        marker=dict(color=color, size=size, opacity=opacity),
        name=label,
    )

def construct_arrow(x, y, z, u=0, v=0, w=1, color='red', size=10, opacity=1.0, label="line"):
    """
    Dibuja una flecha 3D usando plotly.graph_objs.Cone. 
    (u, v, w) es el vector (dirección/magnitude) de la flecha.

    Por defecto apunta en z+.

    Requiere plotly >= 4.x
    """
    return go.Cone(
        x=[x], y=[y], z=[z],
        u=[u], v=[v], w=[w],
        anchor="tail",
        showscale=False,
        sizemode="absolute",
        colorscale=[[0, color], [1, color]],
        opacity=opacity,
        sizeref=size/3,
        name=label,
    )
