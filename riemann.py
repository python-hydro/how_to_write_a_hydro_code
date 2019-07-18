r"""Solve riemann shock tube problem for a general equation of state
using the method of Colella, Glaz, and Ferguson (this is the main
solver used in Castro).  Use a two shock approximation, and linearly
interpolation between the head and tail of a rarefaction to treat
rarefactions.

The Riemann problem for the Euler's equation produces 4 states,
separated by the three characteristics (u - cs, u, u + cs):


        l_1      t    l_2       l_3
         \       ^   .       /
          \  *L  |   . *R   /
           \     |  .     /
            \    |  .    /
        L    \   | .   /    R
              \  | .  /
               \ |. /
                \|./
       ----------+----------------> x

       l_1 = u - cs   eigenvalue
       l_2 = u        eigenvalue (contact)
       l_3 = u + cs   eigenvalue

       only density jumps across l_2

  References:

   CG:   Colella & Glaz 1985, JCP, 59, 264.

   CW:   Colella & Woodward 1984, JCP, 54, 174.

   Fry:  Fryxell et al. 2000, ApJS, 131, 273.

   Toro: Toro 1999, ``Riemann Solvers and Numerical Methods for Fluid
         Dynamcs: A Practical Introduction, 2nd Ed.'', Springer-Verlag

"""

import numpy as np

URHO = 0
UMX = 1
UENER = 2

QRHO = 0
QU = 1
QP = 2

NVAR = 3


def riemann(q_l, q_r, gamma):
    """solve the Riemann problem given left and right primitive variable
    states.  We return the flux"""

    flux = np.zeros(NVAR)

    small_rho = 1.e-10
    small_p = 1.e-10
    small_u = 1.e-10

    rho_l = max(q_l[QRHO], small_rho)
    u_l = q_l[QU]
    p_l = max(q_l[QP], small_p)

    rho_r = max(q_r[QRHO], small_rho)
    u_r = q_r[QU]
    p_r = max(q_r[QP], small_p)

    # wave speeds (Lagrangian sound speed)
    w_l = np.sqrt(gamma*p_l*rho_l)
    w_r = np.sqrt(gamma*p_r*rho_r)

    # construct our guess at pstar and ustar
    wwinv = 1.0/(w_l + w_r)
    p_star = ((w_r*p_l + w_l*p_r) + w_l*w_r*(u_l - u_r))*wwinv
    u_star = ((w_l*u_l + w_r*u_r) + (p_l - p_r))*wwinv

    p_star = max(p_star, small_p)

    if u_star > 0:
        rho_o = rho_l
        u_o = u_l
        p_o = p_l

    elif u_star < 0:
        rho_o = rho_r
        u_o = u_r
        p_o = p_r

    else:
        rho_o = 0.5*(rho_l + rho_r)
        u_o = 0.5*(u_l + u_r)
        p_o = 0.5*(p_l + p_r)

    rho_o = max(rho_o, small_rho)

    c_o = np.sqrt(gamma*p_o/rho_o)

    # compute the rest of the star state
    drho = (p_star - p_o)/c_o**2
    rho_star = rho_o + drho

    c_star = np.sqrt(gamma * p_star/rho_star)

    # sample the solution
    sgn = np.sign(u_star)
    spout = c_o - sgn*u_o
    spin = c_star - sgn*u_star

    ushock = 0.5*(spin + spout)

    if p_star > p_o:
        # compression -- we are a shock
        spin = ushock
        spout = ushock

    if spout-spin == 0.0:
        cavg = 0.5 * (np.sqrt(gamma*p_l/rho_l) + np.sqrt(gamma*p_r/rho_r))
        scr = small_u*0.5*cavg
    else:
        scr = spout - spin

    # interpolate for the case we have a rarefaction
    frac = 0.5*(1.0 + (spout + spin)/scr)
    frac = max(0.0, min(1.0, frac))

    rho_int = frac*rho_star + (1.0 - frac)*rho_o
    u_int = frac*u_star + (1.0 - frac)*u_o
    p_int = frac*p_star + (1.0 - frac)*p_o

    # here we are assuming that the rarefaction spans the interface.  Correct that now
    if spout < 0.0:
        rho_int = rho_o
        u_int = u_o
        p_int = p_o

    if spin >= 0.0:
        rho_int = rho_star
        u_int = u_star
        p_int = p_star

    # now compute the fluxes
    flux[URHO] = rho_int*u_int
    flux[UMX] = rho_int*u_int**2 + p_int
    flux[UENER] = u_int*(p_int/(gamma - 1.0) + 0.5*rho_int*u_int**2 + p_int)

    return flux

if __name__ == "__main__":
    q_l = np.array([1.0, 0.0, 1.0])
    q_r = np.array([0.125, 0.0, 0.1])
    gamma = 1.4

    F = riemann(q_l, q_r, gamma)
    print(F)
