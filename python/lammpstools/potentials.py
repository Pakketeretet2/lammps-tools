"""!
This file contains some pair potentials.

\ingroup lammpstools
"""


import lammpstools
import dumpreader
import numpy as np
import math
import sys


# Curryable pair potentials:
class pair_lj:
    """ !Pair potential for LJ. """
    def __init__(self, epsilon, sigma ):
        self.epsilon = epsilon
        self.sigma = sigma

    def energy(self,r):
        "energy"
        c1 = 4.0*self.epsilon
        r2 = r*r
        sr2 = self.sigma*self.sigma/r2
        sr6 = sr2*sr2*sr2
        return c1*sr6*(sr6-1.0)

    def force(self,r):
        "force"
        c1 = 4.0*self.epsilon
        r2 = r*r
        sr2 = self.sigma*self.sigma/r2
        sr6 = sr2*sr2*sr2
        # The force is 4 epsilon*( -12*s^12/r^13 + 6*s^6/r^7 )
        return c1 * sr6*( 12.0*sr6/r - 6.0/r )
    
    def force_prime(self,r):
        "Derivative of force"
        c1 = 4.0*self.epsilon
        r2 = r*r
        sr2 = self.sigma*self.sigma/r2
        sr6 = sr2*sr2*sr2
        # The force is 4 epsilon*( 12*13*s^12/r^14 - 7*6*s^6/r^8 )
        return c1 * sr6*( 13.0*12.0*sr6/r2 - 7.0*6.0/r2 )


class pair_morse:
    def __init__(self, d0, alpha, r0 ):
        self.d0    = d0
        self.alpha = alpha
        self.r0    = r0

        
    """ !Pair potential for LJ. """
    def energy(self,r):
        ar  = -self.alpha*(r-self.r0)
        edr = math.exp( ar )
        return self.d0*edr*(edr-2.0)

    def force(self,r):
        ar  = -self.alpha*(r-self.r0)
        edr = math.exp( ar )
        a = self.alpha
        
        return self.d0*edr*( 2.0*a*edr - 2.0*a )

    def force_prime(self,r):
        ar  = -self.alpha*(r-self.r0)
        edr = math.exp( ar )
        a = self.alpha
        return self.d0*edr*( 4.0*a*a*edr - 2.0*a*a)

    
class truncated:
    def __init__(self, pair_pot, rc, shifted = False):
        self.pair_pot = pair_pot
        self.rc = rc
        self.shifted = shifted
        self.force_prime_dummy = None
        
        if self.shifted:
            self.E_at_rc = self.pair_pot.energy(rc)
        else:
            self.E_at_rc = 0

        if hasattr(pair_pot, "force_prime"):
            self.force_prime_dummy = force_prime_impl

    def energy(self,r):
        if r < self.rc:
            E = self.pair_pot.energy(r)
            if shifted: E -= self.E_at_rc
            return E
        else:
            return 0.0

    def force(self,r):
        if r > self.rc:
            return self.pair_pot.force(r)
        else:
            return 0
    def force_prime(self,r):
        return force_prime_dummy(self,r)
        

    def force_prime_impl(self,r):
        if r < self.rc:
            return pair_pot.force_prime(r)
        else:
            return 0

    
class smooth_linear:
    def __init__(self, pair_pot, rc ):
        self.pair_pot = pair_pot
        self.rc       = rc

        self.E_at_rc  = pair_pot.energy(rc)
        self.f_at_rc  = pair_pot.force(rc)

    def energy(self,r):
        if r < self.rc:
            return self.pair_pot.energy(r) - self.E_at_rc + self.f_at_rc*(r-self.rc)
        else:
            return 0

    def force(self,r):
        if r < self.rc:
            return self.pair_pot.force(r) - self.f_at_rc
        else:
            return 0

    def force_prime(self,r):
        return self.pair_pot.force_prime(r)



def compute_potential_energy( pair_pot, b ):
    """ !Computes total potential energy and per-atom potential energy
         of entire block. Returns b.meta.th the total and the per-atom energies. """

    e_per_atom = np.zeros(b.meta.N)
    e_total = 0.0
    for i in range(0,b.meta.N):
        for j in range(0,i):
            dist, r = lammpstools.distance( b, i, j )
            e_ij = pair_pot.energy( dist )
            e_per_atom[i] += 0.5*e_ij
            e_per_atom[j] += 0.5*e_ij

            e_total += e_ij
    return e_total, e_per_atom


def test_potential( pair_pot, xlo, xhi, h, tol = 1e-8 ):
    """ !Tests if given potential is consistent using finite differences. """

    problem = 0
    # Test energy/force:
    x = xlo
    while x <= xhi:
        xp = x + h
        xm = x - h
        ep = pair_pot.energy(xp)
        em = pair_pot.energy(xm)
        ff = -(ep - em)/(2.0*h)
        f  = pair_pot.force(x)

        df = ff - f
        fs = ff + f
        
        if( fs < 1e-4 ):
            rel_err = math.fabs( df )
        else:
            rel_err = math.fabs( df ) / math.fabs( 0.5*(fs) )

        if math.fabs( rel_err ) > tol:
            print >> sys.stderr, "Difference f too large! |%e| > %e at %f!" % (rel_err, tol, x)
            if( problem > -1 ): problem -= 1
        x += h

    # Test force/force_prime:
    if not hasattr(pair_pot, "force_prime"):
        return problem

    
    x = xlo
    while x <= xhi:
        xp = x + h
        xm = x - h
        fp = pair_pot.force(xp)
        fm = pair_pot.force(xm)
        ff = -(fp - fm)/(2.0*h)
        f  = pair_pot.force_prime(x)

        df = ff - f
        fs = ff + f
        
        if( fs < 1e-4 ):
            rel_err = math.fabs( df )
        else:
            rel_err = math.fabs( df ) / math.fabs( 0.5*(fs) )
        
        if math.fabs( rel_err ) > tol:        
            print >> sys.stderr, "Difference f' too large! |%e| > %e at %f!" % (df, tol, x)
            if( problem > -2 ): problem -= 2
        x += h

    return problem
