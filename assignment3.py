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
    return t**2


def hermiteCubicSolve(aX, bX, aY, bY, aPrime, bPrime, x):

    term1 = (1+ 2*((x-aX)/(bX - aX))) * ((bX-x)/(bX-aX))**2 * aY
    term2 = (x-aX) * ((bX-x)/(bX-aX))**2 * aPrime
    term3 = (1+ 2*((bX-x)/(bX - aX))) * ((aX-x)/(aX-bX))**2 * bY
    term4 = (x-bX) * ((bX-x)/(aX-bX))**2 * bPrime

    return (term1 + term2 + term3 + term4)


class ComputeYnew:
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

class ComputeError:
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


# a = ComputeYnew(6, 2, 3)
#
# print(a.Ynew)
#
# print(a.k1)
#
# b =ComputeError(0, 0.5, 1, 2, 0.5)
#
# print(b.testMethod())



timeFrame = [0,2]
time = timeFrame[0]
endTime = timeFrame[1]
Ynow = 1
dt = 0.1
tolerance = 0.01
dtmin = 0.05
agrow = 1.25
ashrink = 0.8

tol = 0.01



while time < endTime:
    tryStep = ComputeYnew(time, Ynow, dt)
    errorStep = ComputeError(time, dt, Ynow, tryStep.Ynew, tryStep.k1)
    ei = abs(tryStep.dY - errorStep.dY_alt)/ dt


    if (ei < tol and dt > dtmin):
        time += dt
        Ynow = tryStep.Ynew

        '''altering dt for better estimations'''
        if ei < tol/4:
            dt *= agrow
        elif ei > 0.75* tol:
            dt *= ashrink

        '''making sure dt fits at the end of tolerance'''
        if time+dt > endTime:
            dt = endTime-time;
        elif time + 2*dt > endTime:
            dt = (endTime-time)/2

    else:
        dt /= 2

    print('newStep')
    print(ei)
    print(tryStep.Ynew)
    print(errorStep.dY_alt)



class simTime:

    time = 0.0
    dt = 0.0
    tol = 0.0
    agrow = 0.0
    ashrink = 0.0
    dtmin = 0.0
    dtmax = 0.0
    endTime = 0.0

    stepsSinceRejection = 0
    stepsRejectbed = 0
    stepsAccepted = 0

    """docstring for simTime."""
    def __init__(self, tolerance):

        self.time = 0.0
        self.dt = 0.0
        self.tol = tolerance
        self.agrow = 0.0
        self.ashrink = 0.0
        self.dtmin = 0.0
        self.dtmax = 0.0
        self.endTime = 0.0

        self.stepsSinceRejection = 0
        self.stepsRejectbed = 0
        self.stepsAccepted = 0





#
# def solve(a, b, Ynow, dt):
#
# class Solve(object):
#
#
#     def __init__(self, a, b, Ynow, dt):
#         t_now
#




    # Finally, you should construct an error indicator
    #
    #   ei = || dY - dY_alt ||/dt,
