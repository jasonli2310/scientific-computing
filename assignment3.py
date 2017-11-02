'''
Author: Jason Li
Assignment: CMSC 285 Project 3 - Runge Kutta 4
Professor: DuPont
Date: 10/31/2017
'''

import matplotlib.pyplot as plt
import numpy as np
import math




# You are to implement a solver that approximates the solution
# of ordinary differential systems using a forth order Runge
# Kutta technique. For the problem on [a,b]
#
#    y' = f(t,y),
#    y(a) = y0,


# Given: a, b, f, y0, t0 = a, tn = b,
#
#
# Taking A Step. To advance from knowing the solution up to
# t_now to knowing it up to t_new you will compute four
# derivative values k1, k2, k3, and k4 as follows:
#
#    k1 = f(t_now,Ynow)
#    k2 = f(t_now+dt/2, Ynow + k1 dt/2)
#    k3 = f(t_now+dt/2, Ynow + k2 dt/2)
#    k4 = f(t_now+dt, Ynow + k3 dt)
#    dY = (k1 + 2 k2 + 2 k3 + k4 ) dt/6
#    Ynew = Ynow + dY


def f(t, y):
    return t*y


def hermiteCubicSolve(aX, bX, aY, bY, aPrime, bPrime, x):

    term1 = (1+ 2*((x-aX)/(bX - aX))) * ((bX-x)/(bX-aX))**2 * aY
    term2 = (x-aX) * ((bX-x)/(bX-aX))**2 * aPrime
    term3 = (1+ 2*((bX-x)/(bX - aX))) * ((aX-x)/(aX-bX))**2 * bY
    term4 = (x-bX) * ((bX-x)/(aX-bX))**2 * bPrime

    return (term1 + term2 + term3 + term4)


class computeYnew:
    """docstring for step."""
    def __init__(self, t_now, Ynow, dt):
        self.k1 = f(t_now, Ynow)
        self.k2 = f(t_now + dt/2, Ynow + self.k1 * dt/2)
        self.k3 = f(t_now + dt/2, Ynow + self.k2 * dt/2)
        self.k4 = f(t_now + dt, Ynow + self.k3 * dt)
        self.dY = (self.k1 + 2*self.k2+ 2*self.k2 + self.k4) * dt/2
        self.Ynew = Ynow + self.dY

# You should construct an error indicator in the following
# way.  Once a guess at Ynew exists evaluate f at (t_now+dt,
# Ynew), call this result f_new.  View Y as an Hermite cubic
# on the interval [t_now,t_now+dt] using k1 and f_new as the
# derivatives of Y at the ends, together with the values Ynow
# and Ynew. Evaluate Y at two points

class computeError:
    """docstring for Error."""
    def __init__(self, t_now, dt, Ynow, Ynew, k1):

        self.t_new = t_now+dt

        self.f_new = f(self.t_new,Ynew)
        self.c1 = 0.4-math.sqrt(0.06)
        self.c2 = 0.4+math.sqrt(0.06)

        self.t_c1 = t_now + dt*self.c1
        self.t_c2 = t_now + dt*self.c2

        self.Yc1 = hermiteCubicSolve(t_now, self.t_new, Ynow, Ynew, k1, self.f_new, self.t_c1)
        self.Yc2 = hermiteCubicSolve(t_now, self.t_new, Ynow, Ynew, k1, self.f_new, self.t_c2)

        self.Fc1 = f(self.t_c1, self.Yc1)
        self.Fc2 = f(self.t_c2, self.Yc2)

        self.w1 = (3*self.c2-1)/(6*(self.c1-self.c2)*(self.c1-1))
        self.w2 = (3*self.c1-1)/(6*(self.c2-self.c1)*(self.c2-1))
        self.wnew = 1 - self.w1 - self.w2

        self.dY_alt = dt * ( self.w1 * self.Fc1 + self.w2 * self.Fc2 + self.wnew * self.f_new)


a = computeYnew(6, 2, 3)


# print(f(a.start, a.end))

print(a.Ynew)

print(a.k1)


b =computeError(0, 0.5, 1, 2, 0.5)
print(b.dY_alt)




    # where

    #
    # Call these values Yr1 and Yc2.  Form FM1 =f(t_c1,Yc1), and
    # fc2 =f(t_c2,Yc2). The alternate value of
    #
    # dY_alt = dt( w1 fc1 + w2 fc2 + wnew f_new)
    #
    # is formed by integrating f on the interval (t_now,tnew) using fc1,
    # fc2, and f_new. This integration should be correct on quartics.  With
    # some positive probability, the formulas for the w1, w2, and wnew are
    #
    #     w1 = (3c2-1)/(6(c1-c2)(c1-1)),
    #     w2 = (3c1-1)/(6(c2-c1)(c2-1)),
    #     wnew = 1 - w1 -w2.
    #
    # Finally, you should construct an error indicator
    #
    #   ei = || dY - dY_alt ||/dt,
