import meshio
import numpy as np
from numpy import sin, cos, pi
import plotly.graph_objects as go
from typing import List, Dict

def model_to_figure_mesh3d(
    model_dict: Dict[str, str], # obj/ply file path
    paper_bgcolor: str = '#fff',
    mesh_color: str = '#ddd',
    flatshading: bool = False,
    width: int = 700,
    height: int = 700
):
    
    legend_i=2
    data_traces = []
    for name, path in model_dict.items():
        mesh = meshio.obj.read(path)
        vertices = mesh.points
        triangles = mesh.cells_dict['triangle']

        # The 3d triangulated object, representing Beethoven bust, is represented in a system of 
        # coordinates different from the world system, and is displayed as lying down. 
        # That's why we apply some transformations, such that to see its face with the default settings 
        # of the Plotly camera posit

        def rot_x(t):
            return np.array([[1, 0, 0],
                            [0, cos(t), -sin(t)],
                            [0, sin(t), cos(t)]])
        def rot_z(t):
            return np.array([[cos(t), -sin(t), 0],
                            [sin(t), cos(t), 0],
                            [0, 0, 1]])

        #apply rot_z(angle2) * rot_x(angle1)
        A = rot_x(pi/4)
        B = rot_z(4*pi/9+pi/4)

        #Apply the product of the two rotations to the object vertices:
        new_vertices = np.einsum('ik, kj -> ij',  vertices, (np.dot(B, A)).T)#new_vertices has the shape (n_vertices, 3)

        x, y, z = new_vertices.T
        I, J, K = triangles.T

        tri_points = new_vertices[triangles] 

        colorscale=[0, mesh_color], [1., mesh_color]

        legend_i += 1
        pl_mesh = go.Mesh3d(
            x=x, y=y, z=z,
            colorscale=colorscale, 
            intensity=z,
            flatshading=flatshading, # 平滑
            i=I, j=J, k=K,
            showscale=False,
            name=name,
            legend=f'legend{legend_i}'
        )

        pl_mesh.update(
            cmin=-7,# atrick to get a nice plot (z.min()=-3.31909)
            lighting=dict(
                ambient=0.7, # 环境光
                diffuse=1,
                fresnel=0.1,
                specular=0,
                roughness=1,
                facenormalsepsilon=1e-1,
                vertexnormalsepsilon=1e-15
            ),
            # lightposition=dict(x=100, y=200,z=0)
        )

        Xe, Ye, Ze = [], [], []
        for T in tri_points:
            Xe.extend([T[k%3][0] for k in range(4)]+[ None])
            Ye.extend([T[k%3][1] for k in range(4)]+[ None])
            Ze.extend([T[k%3][2] for k in range(4)]+[ None])

        #define the trace for triangle sides
        lines = go.Scatter3d(
            x=Xe, y=Ye, z=Ze,
            mode='lines',
            name=name,
            line=dict(color='#555', width=1)
        )
        
        data_traces.extend([pl_mesh, lines])

    layout = go.Layout(
    #  title="Beethoven<br>Mesh3d with flatshading",
        font=dict(size=12),
        width=width,
        height=height,
        scene_xaxis_visible=False,
        scene_yaxis_visible=False,
        scene_zaxis_visible=False,
        paper_bgcolor=paper_bgcolor,
        margin=dict(l=0, r=0, t=0, b=0)
    )

    fig = go.Figure(data=data_traces, layout=layout)

    return fig

if __name__ == '__main__':
    path_dict = {
        'M6_10_E7': '/rad/share/omics-viewer/3D_model/Model6_10_E7.obj',
        'M6_9_E3':  '/rad/share/omics-viewer/3D_model/Model6_9_E3.obj'
    }
    fig = model_to_figure_mesh3d(path_dict)
    fig.show()