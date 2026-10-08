"""Complete (axial + radial) induced velocity of semi-infinite right vortex cylinders
(Branlard & Gaunaa, Wind Energy 18, 2015), double precision.
Convention: z is the axial coordinate, positive DOWNSTREAM; cylinder k of radius Rk starts at z = zk
and extends to z -> +infinity. gamma_t < 0 for a wind turbine (u_z < 0 = deficit)."""
import numpy as np
from scipy.special import ellipk, ellipe, elliprf, elliprj

def ellippi(n, m):
    """complete elliptic integral of the third kind Pi(n|m), n < 1, m < 1 (Carlson forms)"""
    return elliprf(0., 1.-m, 1.) + n/3.*elliprj(0., 1.-m, 1., 1.-n)

def ur_cyl(r, z, R, gam):
    """radial velocity at (r, z), z relative to the cylinder start (symmetric in z)"""
    k2 = 4*r*R/((R+r)**2 + z**2); k = np.sqrt(k2)
    return -gam/(2*np.pi)*np.sqrt(R/r)*((2-k2)/k*ellipk(k2) - 2/k*ellipe(k2))

def uz_cyl(r, z, R, gam):
    """axial velocity at (r, z), z relative to the cylinder start; r != R"""
    r = np.asarray(r, float); z = np.asarray(z, float)
    k2 = 4*r*R/((R+r)**2 + z**2); k = np.sqrt(k2)
    k02 = 4*r*R/(R+r)**2
    inside = np.where(r < R, 1.0, 0.0)
    return gam/2*(inside + z*k/(2*np.pi*np.sqrt(r*R))*(ellipk(k2) + (R-r)/(R+r)*ellippi(k02, k2)))

def sweep(rm, zm, Rk, zk, gam):
    """u_r and u_z at the points (rm, zm) from all cylinders (Rk, zk, gam)"""
    ur = np.zeros_like(rm); uz = np.zeros_like(rm)
    for R, zc, g in zip(Rk, zk, gam):
        if abs(g) < 1e-14 or R < 1e-6: continue
        ur += ur_cyl(rm, zm-zc, R, g); uz += uz_cyl(rm, zm-zc, R, g)
    return ur, uz

def strengths(a_inf, U0):
    """gamma_t(k) = 2 U0 (a_k - a_{k-1}), ghost values 0 (Li et al. 2022, Eq. 20)"""
    ae = np.concatenate([[0.], a_inf, [0.]])
    return 2*U0*(ae[1:]-ae[:-1])

if __name__ == '__main__':
    # verification against direct Biot-Savart integration of vortex rings
    from scipy.integrate import quad
    def ring(r, z, R):   # velocity of a ring of unit circulation at axial offset z (point minus ring)
        k2 = 4*r*R/((R+r)**2+z**2); K=ellipk(k2); E=ellipe(k2); d=np.sqrt((R+r)**2+z**2)
        uz = 1/(2*np.pi*d)*(K + (R*R-r*r-z*z)/((R-r)**2+z**2)*E)
        ur = -z/(2*np.pi*r*d)*(K - (R*R+r*r+z*z)/((R-r)**2+z**2)*E)
        return ur, uz
    R=1.0; g=-0.7
    for r,z in [(0.5,0.3),(0.5,-0.3),(0.9,0.05),(1.3,0.4),(1.3,-0.4),(0.97,-0.02)]:
        # the ring formula has the opposite orientation convention -> compare magnitudes via known limits
        ur = quad(lambda s: ring(r, z-s, R)[0], 0, np.inf, limit=400)[0]*g
        uz = quad(lambda s: ring(r, z-s, R)[1], 0, np.inf, limit=400)[0]*g
        print(f"r={r} z={z:+.2f}  ur quad {ur:+.6f} formula {ur_cyl(r,z,R,g):+.6f} | uz quad {uz:+.6f} formula {float(uz_cyl(r,z,R,g)):+.6f}")
    print('inside z=0:', float(uz_cyl(0.4,0.0,R,g)), ' far wake:', float(uz_cyl(0.4,50.,R,g)), ' upstream:', float(uz_cyl(0.4,-50.,R,g)), 'outside z=0', float(uz_cyl(1.4,0.0,R,g)))
