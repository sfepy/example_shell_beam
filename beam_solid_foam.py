"""
Display geometry and results:

    sfepy-view beam_solid_foam.vtk -f mat_id:e:p0 --color-map=cool --camera-position="-0.18,0.07,0.57,0.04,-0.03,0.24,0.12,0.98,-0.19"
    sfepy-view results/beam_solid_foam.vtk -f u:wu:f10:p0 0:vw:p0 --camera-position="-0.4,0.16,0.57,0.02,0,0.23,0.22,0.97,-0.11"
"""
import numpy as nm
import os.path as osp
from sfepy.base.base import Struct
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.discrete import Integral

wdir = osp.dirname(__file__)

filename_mesh, np = 'beam_solid_foam.vtk', 13


def pp_hook(out, pb, state, extend=False):
    stress = pb.evaluate('ev_cauchy_stress.i.Omega(solid.D, u)', mode='el_eval')
    strain = pb.evaluate('ev_cauchy_strain.i.Omega(u)', mode='el_eval')
    out['E'] = Struct(name='strain', mode='cell', data=strain)
    out['S'] = Struct(name='stress', mode='cell', data=stress)

    return out


regions = {
    'Omega': 'all',
    'Left': ('vertices in (z < 0.0001)', 'vertex', 'Omega'),
    'Right': ('vertices in (z > 0.3999)', 'vertex', 'Omega'),
    'Top': ('vertices in (y > 0.01499)', 'vertex', 'Omega'),
    'Edge': ('r.Right *v r.Top', 'vertex', 'Omega'),
    'Omega1': 'cells of group 1',
    'Omega2': 'cells of group 2',
}

force = 1e3
pload = nm.array([[0.0, -force / np, 0.0]] * np)

options = {
    'output_dir': osp.join(wdir, 'results'),
    'post_process_hook': pp_hook,
}

materials = {
    'solid': ({'D': {
        'Omega1': stiffness_from_youngpoisson(3, 210e9, 0.3),
        'Omega2': stiffness_from_youngpoisson(3, 20e9, 0.25),
    }},),
    'load': ({'.val': pload},)
}

fields = {
    'displacement': ('real', 'vector', 'Omega', 1),
}

integrals = {
    'i': 2,
}

variables = {
    'u': ('unknown field', 'displacement', 0),
    'v': ('test field', 'displacement', 'u'),
}

ebcs = {
    'Fixed': ('Left', {'u.all': 0.0}),
}

equations = {
    'balance_of_forces':
    """dw_lin_elastic.i.Omega(solid.D, v, u)
     = dw_point_load.i.Edge(load.val, v)""",
}

solvers = {
    'ls': ('ls.mumps', {}),
    'newton': ('nls.newton', {'eps_a': 1e-5}),
}
