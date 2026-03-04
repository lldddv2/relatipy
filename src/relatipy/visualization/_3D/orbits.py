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

def construct_basic_path_plot(R_s, path, color_path='red', color_black_hole='black', opacity_path=0.2, opacity_black_hole=0.3):
    max_range = np.max(np.abs(path.state_vector[1:]))
    fig = go.Figure(data=[
        construct_circular_plane_surface(max_range*1.2, color_black_hole, opacity_black_hole/3),
        construct_black_hole(R_s, color_black_hole, opacity_black_hole),
        *construct_orbit_path(path, color_path, opacity_path),
    ])

    set_axis_equal(fig, max_range)
    return fig