"""
Display geometry and results:

    sfepy-view beam_shell.vtk -f mat_id:e:p0 --color-map=cool --camera-position="-0.18,0.07,0.57,0.04,-0.03,0.24,0.12,0.98,-0.19"
    sfepy-view results/beam_shell.vtk -f u_disp:wu_disp:f10:p0 0:vw:p0 --camera-position="-0.4,0.16,0.57,0.02,0,0.23,0.22,0.97,-0.11"
"""
import numpy as nm
import os.path as osp
from sfepy.base.base import Struct
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.discrete import Integral
import sfepy.mechanics.shell10x as sh

wdir = osp.dirname(__file__)

filename_mesh, np = 'beam_shell.vtk', 13
thickness = 0.002


def pp_hook(out, pb, state, extend=False):
    return out


regions = {
    'Omega': 'cells of group 1',
    'Left': ('vertices in (z < 0.0001)', 'vertex'),
    'Right': ('vertices in (z > 0.3999)', 'vertex'),
    'Top': ('vertices in (y > 0.01399)', 'vertex'),
    'Edge': ('r.Right *v r.Top', 'vertex'),
}

force = 1e3
pload = nm.array([[0.0, -force / np, 0.0, 0.0, 0.0, 0.0]] * np)

options = {
    'output_dir': osp.join(wdir, 'results'),
    'post_process_hook': pp_hook,
}

materials = {
    'solid': ({
        'D': sh.create_elastic_tensor(young=210e9, poisson=0.3),
        '.drill': 1e-7,
    },),
    'load': ({'.val': pload},)
}

fields = {
    'fu': ('real', 6, 'Omega', 1, 'H1', 'shell10x'),
}

variables = {
    'u': ('unknown field', 'fu', 0),
    'v': ('test field', 'fu', 'u'),
}

# Custom integral.
aux = Integral('i', order=3)
qp_coors, qp_weights = aux.get_qp('3_8')
qp_coors[:, 2] = thickness * (qp_coors[:, 2] - 0.5)
qp_weights *= thickness

integrals = {
    'i': ('custom', qp_coors, qp_weights),
}

ebcs = {
    'Fixed': ('Left', {'u.all': 0.0}),
}

equations = {
    'balance_of_forces: shell':
    """dw_shell10x.i.Omega(solid.D, solid.drill, v, u)
     = dw_point_load.i.Edge(load.val, v)""",
}

solvers = {
    'ls': ('ls.mumps', {}),
    'newton': ('nls.newton', {'eps_a': 1e-6}),
}
