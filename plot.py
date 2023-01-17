import numpy as nm
import meshio
import matplotlib.pyplot as plt
from gen_mesh import t

tags = [
    ('solid', 'u'),
    ('shell', 'u_disp'),
    ('solid_foam', 'u'),
    ('shell_foam', 'us_disp'),
]

fig, ax = plt.subplots()
ax.set_xlabel('z [mm]')
ax.set_ylabel('u_y [mm]')

for tag, var in tags:
    fn = f'results/beam_{tag}.vtk'
    print(fn)
    m = meshio.read(fn)
    pts = m.points
    maxc = nm.max(pts, axis=0)
    mx = maxc[1] if 'shell' in tag else maxc[1] - t/2
    eps = 1e-4
    idxs = nm.where(nm.logical_and(
        nm.logical_and(pts[:, 1] < mx + eps, pts[:, 1] > mx - eps),
        pts[:, 0] == 0))
    z = pts[idxs, 2].squeeze()
    val = m.point_data[var][idxs, 1].squeeze()
    ax.plot(z * 1e3, val * 1e3, label=tag)

ax.grid(True)
ax.legend()

plt.savefig('results.png')
