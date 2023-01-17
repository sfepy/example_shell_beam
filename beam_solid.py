"""
Display geometry and results:

    sfepy-view beam_solid.vtk -f mat_id:e:p0 --color-map=cool --camera-position="-0.18,0.07,0.57,0.04,-0.03,0.24,0.12,0.98,-0.19"
    sfepy-view results/beam_solid.vtk -f u:wu:f10:p0 0:vw:p0 --camera-position="-0.4,0.16,0.57,0.02,0,0.23,0.22,0.97,-0.11"
"""
import numpy as nm
import os.path as osp
from sfepy.base.base import Struct
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.discrete import Integral

wdir = osp.dirname(__file__)

filename_mesh, np = 'beam_solid.vtk', 13


def pp_hook(out, pb, state, extend=False):
    return out


regions = {
    'Omega': 'all',
    'Left': ('vertices in (z < 0.0001)', 'vertex', 'Omega'),
    'Right': ('vertices in (z > 0.3999)', 'vertex', 'Omega'),
    'Top': ('vertices in (y > 0.01499)', 'vertex', 'Omega'),
    'Edge': ('r.Right *v r.Top', 'vertex', 'Omega'),
}

force = 1e3
pload = nm.array([[0.0, -force / np, 0.0]] * np)

options = {
    'output_dir': osp.join(wdir, 'results'),
    'post_process_hook': pp_hook,
}

materials = {
    'solid': ({'D': stiffness_from_youngpoisson(3, young=210e9, poisson=0.3),},),
    'load': ({'.val': pload},),
}

fields = {
    'displacement': ('real', 'vector', 'Omega', 1),
}

variables = {
    'u': ('unknown field', 'displacement', 1),
    'v': ('test field', 'displacement', 'u'),
}

integrals = {
    'is': 2,
}

ebcs = {
    'Fixed': ('Left', {'u.all': 0.0}),
}

equations = {
    'balance_of_forces: solid':
    """dw_lin_elastic.is.Omega(solid.D, v, u)
     = dw_point_load.is.Edge(load.val, v)""",
}

solvers = {
    'ls': ('ls.mumps', {}),
    'newton': ('nls.newton', {'eps_a': 1e-6}),
}
