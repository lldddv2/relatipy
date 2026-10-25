import numpy as np, time, diskdyn as dd, disktab
# (1) ring formula vs brute-force point masses
def brute(R, z, b, H, N=200000):
    ph = (np.arange(N) + 0.5) * 2 * np.pi / N
    dx = R - b * np.cos(ph); dy = -b * np.sin(ph); dz = z
    d3 = (dx * dx + dy * dy + dz * dz + H * H) ** 1.5
    return -(dx / d3).mean(), -(dz / d3).mean()
for (R, z, b, H) in [(0.3, 0.2, 1.0, 0.01), (1.0, 0.003, 1.0, 0.01), (2.5, -0.7, 1.0, 0.0), (1e-4, 0.5, 1.0, 0.0), (0.999, 0.0, 1.0, 0.01)]:
    g = dd.ring_g(R, z, b, H); gb = brute(R, z, b, H)
    print(f"ring R={R} z={z} H={H}: rel err gR {abs(g[0]-gb[0])/max(abs(gb[0]),1e-30):.1e} gz {abs(g[1]-gb[1])/max(abs(gb[1]),1e-30):.1e}")
# (2) quadrature convergence
x20, w20 = np.polynomial.legendre.leggauss(20); x40, w40 = np.polynomial.legendre.leggauss(40)
for (R, z) in [(0.5, 0.0), (0.5, 0.004), (1.0, 0.0), (1.0, 0.01), (0.01, 0.0001), (2.0, 0.1), (0.2, 0.5)]:
    a = dd.disk_g_direct(R, z, 0.0, 1.0, -1.0, 0.01, x20, w20); b = dd.disk_g_direct(R, z, 0.0, 1.0, -1.0, 0.01, x40, w40)
    print(f"quad R={R} z={z}: gR={a[0]:+.6e} gz={a[1]:+.6e}  d20-40: {abs(a[0]-b[0])/abs(b[0]):.1e} {abs(a[1]-b[1])/max(abs(b[1]),1e-30):.1e}")
t0 = time.time(); TR, TZ, tp = disktab.make_table(); print("table", TR.shape, f"{time.time()-t0:.1f}s")
