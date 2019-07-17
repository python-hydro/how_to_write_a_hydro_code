"""Solve riemann shock tube problem for a general equation of state
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
import sys

URHO = 0
UMX = 1
UENER = 2

QRHO = 0
QU = 1
QP = 2

NVAR = 3


def riemann(q_l, q_r, gamma):

    flux = np.zeros(NVAR)

    small_rho = 1.e-10
    small_p = 1.e-10
    small_u = 1.e-10

    rho_l = q_l[QRHO]
    u_l = q_l[QU]
    p_l = max(q_l[QP], smallp)

    rho_r = q_r[QRHO]
    u_r = q_r[QU]
    p_r = max(q_r[QP], smallp)

    # specific volume
    tau_l = 1./max(rho_l, smlrho)
    tau_r = 1./max(rho_r, smlrho)

    # wave speeds (Lagrangian sound speed)
    w_l = np.sqrt(gamma*p_l*rho_l)
    w_r = np.sqrt(gamma*p_r*rho_r)

    # construct our guess at pstar and ustar
    wwinv = 1.0/(w_l + w_r)
    pstar = ((w_r*p_l + w_l*p_r) + w_l*w_r*(u_l - u_r))*wwinv
    ustar = ((w_l*u_l + w_r*u_r) + (p_l - p_r))*wwinv

    pstar = max(pstar, small_p)

    if ustar > 0:
        rho_o = rho_l
        u_o = u_l
        p_o = p_l

    elif ustar < 0:
        rho_o = rho_r
        u_o = u_r
        p_o = p_r

    else:
        rho_o = 0.5*(rho_l + rho_r)
        u_o = 0.5*(u_l + u_r)
        p_o = 0.5*(p_l + p_r)

    rho_o = max(rho_o, small_rho)

    co = np.sqrt(gamma*p_o/rho_o)

    # compute the rest of the star state
    drho = (pstar - p_o)/co**2
    rho_star = rho_o + drho

    cstar = np.sqrt(gamma * pstar/rho_star)

    # sample the solution
    sgn = np.sign(ustar)
    spout = co - sgn*uo
    spin = cstar - sgn*ustar

    ushock = 0.5*(spin + spout)

    if pstar > p_o:
        # compression -- we are a shock
        spin = ushock
        spout = ushock

    if 

    # now compute the fluxes
    flux[URHO] = rhoav*uav
    flux[UMX] = rhoav*uav*uav + pav
    flux[UENER] = uav*(pav/(gamma - 1.0) + 0.5*rhoav*uav*uav + pav)

    return flux
