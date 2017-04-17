"""! 
Simple quaternion manipulation.

\ingroup lammpstools
"""

import numpy as np
import math
import sys
import copy

class quat:
    """ A simple unit quaternion. """

    def __init__(self,r = 1.0, u = 0.0, v = 0.0, w = 0.0):
        self.r = r
        self.u = u
        self.v = v
        self.w = w

    def __add__(self, other):

        qret = quat( self.r + other.r,
                     self.u + other.u,
                     self.v + other.v,
                     self.w + other.w )
        
        return qret

    def __sub__(self,other):

        qret = quat( self.r - other.r,
                     self.u - other.u,
                     self.v - other.v,
                     self.w - other.w )
        
        return qret

    def __mul__( self, other ):
        # tougher, quaternion multiplication:
        a1 = copy.deepcopy(self.r)
        b1 = copy.deepcopy(self.u)
        c1 = copy.deepcopy(self.v)
        d1 = copy.deepcopy(self.w)
        a2 = copy.deepcopy(other.r)
        b2 = copy.deepcopy(other.u)
        c2 = copy.deepcopy(other.v)
        d2 = copy.deepcopy(other.w)

        qret = quat(
            a1*a2 - b1*b2 - c1*c2 - d1*d2,
            a1*b2 + b1*a2 + c1*d2 - d1*c2,
            a1*c2 - b1*d2 + c1*a2 + d1*b2,
            a1*d2 + b1*c2 - c1*b2 + d1*a2
        )

        return qret

    def __div__( self, scalar ):
        # tougher, quaternion multiplication:
        self.r /= scalar
        self.u /= scalar
        self.v /= scalar
        self.w /= scalar

        
        return self


    
    def __getitem__( self, idx ):
        if idx < 0 or idx > 3:
            raise RuntimeError("Index out of bounds!")

        if idx == 0:
            return self.r
        elif idx == 1:
            return self.u
        elif idx == 2:
            return self.v
        elif idx == 3:
            return self.w
        


    def t(self):
        """ Complex conjugate """
        qret = quat( self.r,
                     -self.u,
                     -self.v,
                     -self.w )
        return qret

    def norm(self):
        r = self.r
        u = self.u
        v = self.v
        w = self.w
        return r*r + u*u + v*v + w*w

    

    def print_me(self):
        print("[ ", self.r, " ", self.u, " ", self.v, " ", self.w, " ]")


    def init_from_spherical(self, x, y, z):
        R2 = x*x + y*y + z*z
        R  = math.sqrt(R2)
        self.r = 0
        self.u = x / R
        self.v = y / R
        self.w = z / R

        return self



def rot_vector( v, axis, angle ):
    """ Rotates given vector v about axis with given angle. """

    # Normalize the vectors:
    axis /= np.linalg.norm(axis)

    q1 = quat( 0, v[0], v[1], v[2] )
    
    s = math.sin(angle/2.0)
    c = math.cos(angle/2.0)
    
    rotator = quat( c, axis[0]*s, axis[1]*s, axis[2]*s )
   
    qrr = rotator * q1
    rotator = rotator.t()
    qr = qrr*rotator

    # Check:
    if abs(qr[0]) > 1e-8:
        print("Real part of rotated quaternion too large! ", qr[0],
              " > ", 1e-8, file = sys.stderr)
    
    vrot = np.array( [ qr[1], qr[2], qr[3] ] )

    return vrot
    
        
    
def test_quat():
    """ !Simple test for quaternions. """

    q1 = quat()
    q2 = quat(1,2,3,4)

    q3 = q1 + q2
    q4 = q1 - q2

    q1.print_me()
    q2.print_me()
    q3.print_me()
    q4.print_me()
    
    q5 = quat().init_from_spherical( 0, 0, 1 )
    q6 = quat().init_from_spherical( 0, 1, 0 )
    q5.print_me()
    q6.print_me()

    # To see if it works as expected, rotate q5 over q6 with
    # 90 degrees, should get -e_x.
    print("(0, 0, 1) about (0, 1, 0) with pi")
    v = rot_vector( np.array( [ 0., 0, 1 ] ), np.array( [ 0., 1, 0 ] ), math.pi )
    print("Got       %f %f %f" % (v[0], v[1], v[2]))
    print("Should be %f %f %f" % (0, 0, -1))

    print("")
    print("(0, 0, 1) about (0, 1, 0) with pi/2")
    v = rot_vector( np.array( [ 0., 0, 1 ] ), np.array( [ 0., 1, 0 ] ), math.pi/2.0 )
    print("Got       %f %f %f" % (v[0], v[1], v[2]))
    print("Should be %f %f %f" % (1,0,0))

    print("")
    print("(0, 0, 1) about (1, 1, 0) with pi")
    v = rot_vector( np.array( [ 0, 0, 1. ] ), np.array( [ 1., 1, 0 ] ), math.pi )
    print("Got       %f %f %f" % (v[0], v[1], v[2]))
    print("Should be %f %f %f" % (0, 0, -1))
    
