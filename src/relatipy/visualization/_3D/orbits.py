from .base_elements import construct_path, construct_point, construct_arrow, construct_sphere_surface, construct_circular_plane_surface
import plotly.graph_objects as go
import numpy as np

def construct_orbit_path(path, color='red', opacity=1.0, size=5, start_point=True, end_point=True, label="orbit"):
    path = path.convert_to("Cartesian")

    elements = [
        construct_path(path, color, opacity, label=label),
    ]
    if start_point:
        start_opacity = opacity + 0.2 if opacity < 0.8 else 1.0
        elements.append(
            construct_point(
                path[1, 0], 
                path[2, 0], 
                path[3, 0], 
                color=color, 
                size=size, 
                opacity=start_opacity,
                label=f"{label} start point",
            ),
        )
    if end_point:
        end_opacity = opacity + 0.2 if opacity < 0.8 else 1.0
        elements.append(
            construct_arrow(
                path[1, -1], 
                path[2, -1],
                path[3, -1],
                u=path[4, -1],
                v=path[5, -1],
                w=path[6, -1],
                color=color,
                size=size,
                opacity=end_opacity,
                label=f"{label} end point",
            ),
        )

    return elements

def construct_black_hole(R_s, color='black', opacity=1.0):
    return construct_sphere_surface(R_s, color, opacity)

def set_axis_equal(fig, max_range):
    fig.update_layout(
        scene=dict(
                xaxis=dict(range=[-max_range, max_range]),
                yaxis=dict(range=[-max_range, max_range]),
                zaxis=dict(range=[-max_range, max_range]),
                aspectmode="cube"
            )
    )

def construct_basic_path_plot(R_s, path, color_path='red', color_black_hole='black', opacity_path=0.2, opacity_black_hole=0.3, plot_plane=True, plot_black_hole=True, plot_center=True):
    path = path.convert_to("Cartesian")
    max_range = np.max(np.abs(path.state_vector[1:]))

    elements = []
    if plot_plane:
        elements.append(construct_circular_plane_surface(max_range*1.2, color_black_hole, opacity_black_hole/3))
    if plot_black_hole:
        elements.append(construct_black_hole(R_s, color_black_hole, opacity_black_hole))
    if plot_center:
        elements.append(construct_point(0, 0, 0, color='black', size=10, opacity=1.0, label="center"))

    elements.extend(construct_orbit_path(path, color_path, opacity_path))

    fig = go.Figure(data=elements)

    set_axis_equal(fig, max_range*1.2)
    return fig


def construct_isco(r: float, prograde: bool, color: str, n_arrows: int = 8, line_width: int = 4, opacity: float = 0.9) -> list:
    """
    Construye el círculo ISCO con flechas de rotación tangenciales centradas sobre el círculo.

    Parameters
    ----------
    r : float
        Radio del ISCO en unidades geométricas.
    prograde : bool
        Si True, sentido antihorario (pro-grado). Si False, horario (retro-grado).
    color : str
        Color del círculo y las flechas.
    n_arrows : int
        Número de flechas distribuidas sobre el círculo.
    line_width : int
        Grosor de línea en píxeles de pantalla.
    opacity : float
        Opacidad del círculo.

    Returns
    -------
    list[go.Scatter3d]
    """
    label = "↺ ISCO pro-grado" if prograde else "↻ ISCO retro-grado"
    sign = 1 if prograde else -1

    phis = np.linspace(0, 2 * np.pi, 120)
    circle = go.Scatter3d(
        x=r * np.cos(phis), y=r * np.sin(phis), z=np.zeros(120),
        mode='lines',
        line=dict(color=color, width=line_width, dash='dashdot'),
        opacity=opacity,
        name=label,
        legendgroup=label,
    )

    length   = r * 0.3
    head_len = 0.2 * length
    spread   = 0.12 * length

    arrow_phis = np.linspace(0, 2 * np.pi, n_arrows, endpoint=False)
    arrows = []
    for phi in arrow_phis:
        u_hat = sign * (-np.sin(phi))
        v_hat = sign * ( np.cos(phi))

        tip_x = r * np.cos(phi)
        tip_y = r * np.sin(phi)

        bx = tip_x - head_len * u_hat
        by = tip_y - head_len * v_hat

        px, py = -v_hat, u_hat  # perpendicular

        arrows.append(go.Scatter3d(
            x=[tip_x, bx + spread * px, None, tip_x, bx - spread * px],
            y=[tip_y, by + spread * py, None, tip_y, by - spread * py],
            z=[0, 0, None, 0, 0],
            mode='lines',
            line=dict(color=color, width=line_width),
            opacity=1,
            hoverinfo='skip',
            showlegend=False,
            legendgroup=label,
        ))

    return [circle] + arrows