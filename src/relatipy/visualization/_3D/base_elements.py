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


def construct_square_plane_surface(semi_side, color='black', opacity=1.0):
    """Build a square plane z=0 using a 2D grid from -semi_side to semi_side."""
    x_s, y_s = np.meshgrid(
        np.linspace(-semi_side, semi_side, 50),
        np.linspace(-semi_side, semi_side, 50),
    )
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


def construct_surface(x, y, z, color='black', opacity=0.5, name="surface"):
    """Build a generic go.Surface from precomputed x, y, z meshes."""
    surface = go.Surface(
        x=x,
        y=y,
        z=z,
        opacity=opacity,
        colorscale=[[0, color], [1, color]],
        showscale=False,
        name=name,
        showlegend=True,
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


def construct_arrow_scatter(x, y, z, u=0, v=0, w=1, color='red', line_width=4, opacity=1.0, label="arrow"):
    """
    Flecha 3D como Scatter3d: asta (línea) + punta en V (dos segmentos). El grosor es en
    unidades de pantalla, así que se ajusta al hacer zoom (no queda fijo como Cone).
    (u, v, w) es el vector dirección/magnitud de la flecha.
    """
    tip_x, tip_y, tip_z = x + u, y + v, z + w
    length = np.sqrt(u * u + v * v + w * w)
    if length < 1e-12:
        line_trace = go.Scatter3d(
            x=[x], y=[y], z=[z],
            mode='lines',
            line=dict(color=color, width=line_width),
            opacity=opacity,
            name=label,
            legendgroup=label,
            showlegend=True,
        )
        return [line_trace]
    # Unit direction and arrowhead (V) in a plane: symmetric V with one perpendicular
    dx, dy, dz = u / length, v / length, w / length
    head_len = 0.2 * length
    spread = 0.12 * length
    back_x = tip_x - head_len * dx
    back_y = tip_y - head_len * dy
    back_z = tip_z - head_len * dz
    # One unit vector perpendicular to (dx, dy, dz); V wings = back ± spread*perp
    if abs(dz) < 0.9:
        ref = np.array([0, 0, 1])
    else:
        ref = np.array([1, 0, 0])
    perp = np.array([dy * ref[2] - dz * ref[1], dz * ref[0] - dx * ref[2], dx * ref[1] - dy * ref[0]])
    n = np.sqrt(perp[0]**2 + perp[1]**2 + perp[2]**2)
    if n > 1e-12:
        perp = perp / n
    else:
        perp = np.array([1, 0, 0])
    left_x = back_x + spread * perp[0]
    left_y = back_y + spread * perp[1]
    left_z = back_z + spread * perp[2]
    right_x = back_x - spread * perp[0]
    right_y = back_y - spread * perp[1]
    right_z = back_z - spread * perp[2]
    shaft = go.Scatter3d(
        x=[x, tip_x],
        y=[y, tip_y],
        z=[z, tip_z],
        mode='lines',
        line=dict(color=color, width=line_width),
        opacity=opacity,
        name=label,
        legendgroup=label,
        showlegend=True,
    )
    head_left = go.Scatter3d(
        x=[tip_x, left_x],
        y=[tip_y, left_y],
        z=[tip_z, left_z],
        mode='lines',
        line=dict(color=color, width=line_width),
        opacity=opacity,
        legendgroup=label,
        showlegend=False,
    )
    head_right = go.Scatter3d(
        x=[tip_x, right_x],
        y=[tip_y, right_y],
        z=[tip_z, right_z],
        mode='lines',
        line=dict(color=color, width=line_width),
        opacity=opacity,
        legendgroup=label,
        showlegend=False,
    )
    return [shaft, head_left, head_right]
