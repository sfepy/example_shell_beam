import sys
import numpy as nm
import meshio

# dimensions
w, h, l, t = 0.04, 0.03, 0.4, 0.002

# elements per dimension
nw, nh, nl, nt, fname, foam = 12, 8, 40, 0, 'beam_shell_foam.vtk', True
# nw, nh, nl, nt, fname, foam = 12, 8, 40, 0, 'beam_shell.vtk', False
# nw, nh, nl, nt, fname, foam = 12, 8, 40, 8, 'beam_solid_foam.vtk', True
# nw, nh, nl, nt, fname, foam = 12, 8, 40, 8, 'beam_solid.vtk', False


def stack_data(data, n=2):
    if data is not None:
        for k in data.keys():
            v = data[k]
            if len(v.shape) > 1:
                data[k] = nm.vstack([v] * n)
            else:
                data[k] = nm.hstack([v] * n)
    return data


def find_master_slave(nodes, tol=1e-9):
    from scipy.spatial import cKDTree

    tr = cKDTree(nodes)
    mtx = tr.sparse_distance_matrix(tr, tol).tocsr()
    nrow = nm.diff(mtx.indptr)
    idxs = nm.where(nrow > 1)[0]

    npairs_max = nm.sum(nrow[idxs] - 1)

    out = nm.empty((npairs_max, 2), dtype=nm.int64)
    idx0 = 0
    for ii in idxs:
        i1, i2 = mtx.indptr[ii], mtx.indptr[ii + 1]
        cols = mtx.indices[i1:i2]
        if cols[cols < ii].shape[0] == 0:
            nc = cols.shape[0]
            if nc == 2:
                out[idx0, :] = cols
                idx0 += 1
            else:
                idx1 = idx0 + nc - 1
                out[idx0:idx1, 0] = cols[0]
                out[idx0:idx1, 1] = cols[1:]
                idx0 = idx1

    return out[:idx0, :]


def mirror_mesh(nodes, elems, pdata=None, cdata=None, c0=0, dim=0):
    remap = [
        {
            '2_4': nm.array([1, 0, 3, 2]),  # quad
            '2_2': nm.array([1, 0]),  # line
        },
        {
            '2_4': nm.array([3, 2, 1, 0]),  # quad
            '2_2': nm.array([1, 0]),  # line
        }
    ]

    nnd, ndim = nodes.shape
    nodes2 = nodes.copy()
    nodes2[:, dim] = 2 * c0 - nodes2[:, dim]
    mnodes = nm.vstack([nodes, nodes2])

    melems = []
    mcdata = None if cdata is None else []
    for k, eg in enumerate(elems):
        et = f'{ndim}_{eg.shape[1]}'
        melems.append(nm.vstack([eg, (eg + nnd)[:, remap[dim][et]]]))
        if cdata is not None:
            mcdata.append(stack_data(cdata[k]))

    return mnodes, melems, stack_data(pdata), mcdata


def merge_nodes(nodes, elems, pdata=None, tol=1e-9):
    ms_tab = find_master_slave(nodes, tol=tol)
    remap = nm.ones((nodes.shape[0],), dtype=nm.int32)
    remap[ms_tab[:, 1]] = -1
    ndidxs = nm.where(remap > 0)[0]
    remap[ndidxs] = nm.arange(len(ndidxs))
    new_nodes = nodes[ndidxs, :]
    remap[ms_tab[:, 1]] = remap[ms_tab[:, 0]]

    if pdata is not None:
        for k in pdata.keys():
            pdata[k] = pdata[k][ndidxs, ...]

    return new_nodes, [remap[eg] for eg in elems], pdata


def extrude_mesh(nodes, elems, pdata, cdata, z, nz):
    remap = {
        2: nm.array([1, 0, 2, 3]),  # line -> quad
        3: nm.array([]),  # tri -> pyramid
        4: nm.array([1, 0, 4, 5, 2, 3, 7, 6]),  # quad -> hexa
    }

    dz = -z / nz
    nnd, ndim = nodes.shape
    new_nodes = nm.zeros((nnd * (nz + 1), 3), dtype=nm.float64)
    new_nodes[:, :2] = nm.tile(nodes, (nz + 1, 1))
    new_nodes[:, 2] = nm.repeat(nm.arange(nz + 1) / nz * z, nnd)

    new_elems = []
    new_cdata = None if cdata is None else []
    for k, eg in enumerate(elems):
        nel, nen = eg.shape
        nen2 = nen * 2
        new_eg = nm.zeros((nel * nz, nen2), dtype=nm.int64)
        new_eg[:, :nen] = nm.tile(eg, (nz, 1))
        new_eg[:, nen:] = nm.tile(eg + nnd, (nz, 1))
        new_eg = new_eg[:, remap[eg.shape[1]]]
        new_eg += nm.repeat(nm.arange(nz) * nnd, nel)[:, None]
        new_elems.append(new_eg)

        if cdata is not None:
            new_cdata.append(stack_data(cdata[k], nz))

    new_pdata = stack_data(pdata, nz + 1)

    return new_nodes, new_elems, new_pdata, new_cdata


def mesh_vol():
    n1, n2 = nw // 2, nh // 2
    n = n1 + n2

    aux = nm.meshgrid(nm.arange(n + 1), nm.arange(nt + 1))
    nodes = nm.array(aux, dtype=nm.float64).reshape((2, -1)).T
    flag = nm.ones(((n + 1) * (nt + 1),), dtype=nm.int64)
    flag[nodes[:, 0] > n1] = 2

    nodes[flag == 1, 0] *= w / nw
    nodes[flag == 1, 1] *= t / nt
    nodes[flag == 1, 1] += h / 2 - t

    aux = nodes[flag == 2, 0].copy()
    nodes[flag == 2, 0] = nodes[flag == 2, 1] * t / nt + w / 2 - t
    nodes[flag == 2, 1] = (n2 - (aux - n1)) * h / nh

    idxs = nm.arange(nt + 1) * (n + 1) + n1
    y = nm.arange(nt + 1) * t / nt + h / 2 - t
    nodes[idxs, 1] = y
    nodes[idxs, 0] = y + (w - h) / 2

    elems = nm.tile([0, 1, n + 2, n + 1], (n * nt, 1))
    elems += (nm.arange(n * nt) + nm.repeat(nm.arange(nt), n))[:, None]

    mat_id = nm.ones((elems.shape[0],), dtype=nm.int64)
    return nodes, elems, mat_id


def mesh_shell():
    n1, n2 = nw // 2, nh // 2
    n = n1 + n2

    nodes = nm.zeros((n + 1, 2), dtype=nm.float64)
    nodes[:, 0] = nm.arange(n + 1)
    flag = nm.ones((n + 1,), dtype=nm.int64)
    flag[nodes[:, 0] > n1] = 2

    nodes[flag == 1, 0] *= w / nw
    nodes[flag == 1, 1] = h/2 - t/2

    aux = nodes[flag == 2, 0].copy()
    nodes[flag == 2, 0] = w/2 - t/2
    nodes[flag == 2, 1] = (n2 - (aux - n1)) * h / nh

    nodes[n1, 0] = w/2 - t/2
    nodes[n1, 1] = h/2 - t/2

    elems = nm.tile([0, 1], (n, 1))
    elems += nm.arange(n)[:, None]

    mat_id = nm.ones((elems.shape[0],), dtype=nm.int64)

    return nodes, elems, mat_id


def mesh_in(t):
    n1, n2 = nw // 2, nh // 2
    aux = nm.meshgrid(nm.arange(n1 + 1) / n1 * w / 2,
                      nm.arange(n2 + 1) / n2 * h / 2)
    nodes = nm.array(aux, dtype=nm.float64).reshape((2, -1)).T
    nodes[nodes[:, 0] > 0.499 * w, 0] = w / 2 - t
    nodes[nodes[:, 1] > 0.499 * h, 1] = h / 2 - t

    elems = nm.tile([0, 1, n1 + 2, n1 + 1], (n1 * n2, 1))
    elems += (nm.arange(n1 * n2)
              + nm.repeat(nm.arange(n2), n1))[:, None]

    mat_id = nm.ones((elems.shape[0],), dtype=nm.int64) + 1

    return nodes, elems, mat_id


def main():

    pdata = {}

    nodes, elems, mat_id = mesh_shell() if nt == 0 else mesh_vol()

    if foam:
        ti = t / 2 if nt == 0 else t
        nodes2, elems2, mat_id2 = mesh_in(ti)
        elems = [elems, elems2 + nodes.shape[0]]
        nodes = nm.vstack([nodes, nodes2])
        nodes, elems, pdata = merge_nodes(nodes, elems, pdata)
        cdata = [{'mat_id': mat_id}, {'mat_id': mat_id2}]
    else:
        elems = [elems]
        cdata = [{'mat_id': mat_id}]

    nodes, elems, pdata, cdata = mirror_mesh(nodes, elems, pdata, cdata, dim=0)
    nodes, elems, pdata, cdata = mirror_mesh(nodes, elems, pdata, cdata, dim=1)
    nodes, elems, pdata = merge_nodes(nodes, elems, pdata)
    nodes, elems, pdata, cdata = extrude_mesh(nodes, elems, pdata, cdata, l, nl)

    cdata_out = {}
    keys = set([k for cdg in cdata for k in cdg.keys()])
    cdata_out = {k: [cdg[k] for cdg in cdata] for k in keys}

    etypes = {
        8: 'hexahedron',
        4: 'quad',
        2: 'line',
    }

    elems_out = []
    for eg in elems:
        elems_out.append((etypes[eg.shape[1]], eg))

    m = meshio.Mesh(nodes, elems_out, cell_data=cdata_out)

    print(f'>>> {fname} <<<')
    m.write(fname, binary=False)


if __name__ == '__main__':
    sys.exit(main())
