"""
Display geometry and results:

    sfepy-view beam_shell_foam.vtk -f mat_id:e:m2:o0:p0 mat_id:e:m1:o1:p0 --color-map=cool --camera-position="-0.18,0.07,0.57,0.04,-0.03,0.24,0.12,0.98,-0.19"
    sfepy-view beam_shell_foam.vtk -f mat_id:e:m2:o1:p0 mat_id:e:m1:o0:p0 --color-map=cool --camera-position="-0.18,0.07,0.57,0.04,-0.03,0.24,0.12,0.98,-0.19"
    sfepy-view results/beam_shell_foam.vtk -f uf:wuf:f10:m2:p0 0:vw:m2:p0 us_disp:wus_disp:f10:m1:p1 --camera-position="-0.4,0.16,0.57,0.02,0,0.23,0.22,0.97,-0.11" --grid-vector1="0, 1.6, 0"
"""
import numpy as nm
import os.path as osp
from sfepy.base.base import Struct
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.discrete import Integral
import sfepy.mechanics.shell10x as sh

wdir = osp.dirname(__file__)

filename_mesh, np = 'beam_shell_foam.vtk', 13
thickness = 0.002


def match_nodes(coors1, coors2):
    """
    Match coordinates `coors1` with `coors2`.
    """
    from sfepy.discrete.fem.mesh import find_map

    if coors1.shape != coors2.shape:
        raise ValueError('incompatible shapes: %s == %s'\
                         % (coors1.shape, coors2.shape))

    i1, i2 = find_map(coors1, coors2, join=False)

    if i1.shape[0] != coors1.shape[0]:
        print(coors1[i1])
        print(coors2[i2])
        ii = nm.setdiff1d(nm.arange(coors1.shape[0]), i1)
        print(coors1[ii])
        print(coors2[ii])
        raise ValueError('cannot match nodes!')

    return i1, i2


def pp_hook(out, pb, state, extend=False):
    stress = pb.evaluate('ev_cauchy_stress.if.OmegaF(foam.D, uf)',
                         mode='el_eval')
    strain = pb.evaluate('ev_cauchy_strain.if.OmegaF(uf)', mode='el_eval')
    out['E'] = Struct(name='strain', mode='cell', data=strain, region='OmegaF')
    out['S'] = Struct(name='stress', mode='cell', data=stress, region='OmegaF')

    idxs = pb.domain.regions['OmegaS'].vertices

    v = out['us_disp'].data
    aux = nm.zeros((pb.domain.shape.n_nod, v.shape[1]), dtype=v.dtype)
    aux[idxs] = v
    out['us_disp'].data = aux

    v = out['us_rot'].data
    aux = nm.zeros((pb.domain.shape.n_nod, v.shape[1]), dtype=v.dtype)
    aux[idxs] = v
    out['us_rot'].data = aux

    return out


regions = {
    # shell
    'OmegaS': 'cells of group 1',
    'LeftS': ('vertices in (z < 0.0001)', 'vertex', 'OmegaS'),
    'RightS': ('vertices in (z > 0.3999)', 'vertex', 'OmegaS'),
    'TopS': ('vertices in (y > 0.01399)', 'vertex', 'OmegaS'),
    'EdgeS': ('r.RightS *v r.TopS', 'vertex', 'OmegaS'),
    # solid foam
    'OmegaF': 'cells of group 2',
    'LeftF': ('vertices in (z < 0.0001)', 'facet', 'OmegaF'),
    'RightF': ('vertices in (z > 0.3999)', 'facet', 'OmegaF'),
    'SurfaceF_': ('vertices of surface', 'facet', 'OmegaF'),
    'SurfaceF__': ('r.SurfaceF_ -s r.LeftF', 'facet', 'OmegaF'),
    'SurfaceF': ('r.SurfaceF__ -s r.RightF', 'facet', 'OmegaF'),
}

force = 1e3
pload = nm.array([[0.0, -force / np, 0.0, 0.0, 0.0, 0.0]] * np)

options = {
    'output_dir': osp.join(wdir, 'results'),
    'post_process_hook': pp_hook,
}

materials = {
    'shell': ({
        'D': sh.create_elastic_tensor(young=210e9, poisson=0.3),
        '.drill': 1e-7,
    },),
    'load': ({'.val': pload},),
    'foam': ({'D': stiffness_from_youngpoisson(3, 20e9, 0.25), },),
}

fields = {
    'fu': ('real', 6, 'OmegaS', 1, 'H1', 'shell10x'),
    'displacement': ('real', 'vector', 'OmegaF', 1),
}

variables = {
    'us': ('unknown field', 'fu', 0),
    'vs': ('test field', 'fu', 'us'),
    'uf': ('unknown field', 'displacement', 1),
    'vf': ('test field', 'displacement', 'uf'),
}

# Custom integral.
aux = Integral('i', order=3)
qp_coors, qp_weights = aux.get_qp('3_8')
qp_coors[:, 2] = thickness * (qp_coors[:, 2] - 0.5)
qp_weights *= thickness

integrals = {
    'is': ('custom', qp_coors, qp_weights),
    'if': 2,
}

ebcs = {
    'FixedS': ('LeftS', {'us.all': 0.0}),
    'FixedF': ('LeftF', {'uf.all': 0.0}),
}

functions = {
    'match_nodes': (match_nodes,),
}

lcbcs = {
    'match_surface': (['OmegaS', 'SurfaceF'], {'us.[0,1,2]': 'uf.[0,1,2]'},
                      'match_nodes', 'match_dofs'),
}

equations = {
    'balance_of_forces - shell':
    """dw_shell10x.is.OmegaS(shell.D, shell.drill, vs, us)
     = dw_point_load.is.EdgeS(load.val, vs)""",
    'balance_of_forces - foam':
    """dw_lin_elastic.if.OmegaF(foam.D, vf, uf) = 0""",
}

solvers = {
    'ls': ('ls.mumps', {}),
    'newton': ('nls.newton', {'eps_a': 1e-6}),
}
