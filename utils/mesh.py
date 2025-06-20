import plotly.graph_objects as go
from typing import List, Dict
from dash_obj_in_3dmesh import geometry_tools

def obj_mtl_to_mesh3d(
    model_names: str | List[str],
    path: str = '/data1/share/omics-viewer/3D_model',
    lighting: Dict[str, float] = {'ambient': 0.6, 'specular': 0.2, 'diffuse': 1, 'roughness': 0.5},
    lightposition: Dict[str, int] = { 'x': 0, 'y': 0, 'z': 100 },
    height=500,
    width=500,
    paper_bgcolor='#fff'
):
    if isinstance(model_names, str):
        model_names = [model_names]
    data = geometry_tools.import_geometry(model_names, path=path)
    layout = go.Layout(
        font=dict(size=12),
        paper_bgcolor=paper_bgcolor,
        margin=dict(l=0, r=0, t=0, b=0),
        width=width,
        height=height,
        scene_xaxis_visible=False,
        scene_yaxis_visible=False,
        scene_zaxis_visible=False,
    )
    fig = go.Figure(
        data = data, layout = layout
    ).update_traces(
        opacity=1,
        lighting = lighting,
        lightposition = lightposition,
        showlegend=False,
    )
    
    return fig