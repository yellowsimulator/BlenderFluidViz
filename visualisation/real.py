#!/usr/bin/env python
"""
2D wave equation solved by finite differences::

  dt, cpu_time = solver(I, V, f, c, Lx, Ly, Nx, Ny, dt, T,
                        user_action=None, version='scalar',
                        stability_safety_factor=1)

Solve the 2D wave equation u_tt +b*u_t= d_x(q_x(x,y)*u_x) + d_y(q(x,y)*u_y) + f(x,y,t) on (0,L)*(0,L) with
du/dn=0 on the boundary and initial condition du/dt=V(x,y) and u(x,y,0)=I(x,y)

Nx and Ny are the total number of mesh cells in the x and y
directions. The mesh points are numbered as (0,0), (1,0), (2,0),
..., (Nx,0), (0,1), (1,1), ..., (Nx, Ny).

dt is the time step. If dt<=0, an optimal time step is used.
T is the stop time for the simulation.

I, V, f are functions: I(x,y), V(x,y), f(x,y,t). V and f
can be specified as None or 0, resulting in V=0 and f=0.

user_action: function of (u, x, y, t, n) called at each time
level (x and y are one-dimensional coordinate vectors).
This function allows the calling code to plot the solution,
compute errors, etc.
"""
#from scitools import *
#from numpy import *
#from matplotlib import *
#from pylab import *
#from mpl_toolkits.mplot3d import axes3d
#import matplotlib.pyplot as plt
#import numpy as np
#from scitools.easyviz import *
import time, sys
from math import exp
#from scitools import *
from numpy import *
from matplotlib import *
from pylab import *
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np


import time, sys
import os
from scitools.std import *
import numpy


def solver(I, V, f, q, Lx, Ly, Nx, Ny, dt, T, b,
           user_action=None, version='scalar'):
        
    if version == 'vectorized':
        advance = advance_vectorized
    elif version == 'scalar':
        advance = advance_scalar

    x = linspace(0, Lx, Nx+1)  # mesh points in x dir
    y = linspace(0, Ly, Ny+1)  # mesh points in y dir
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    xv = x[:,newaxis]          # for vectorized function evaluations
    yv = y[newaxis,:]
    c = 1
    stability_limit = (1/float(c))*(1/sqrt(1/dx**2 + 1/dy**2)) 
    if dt <= 0:                # max time step?
        safety_factor = -dt    # use negative dt as safety factor
        dt = safety_factor*stability_limit
    elif dt > stability_limit:
        print('error:{} exceeds the stability limit {}'.format(dt,stability_limit))
              
    Nt = int(round(T/float(dt)))
    t = linspace(0, Nt*dt, Nt+1)    # mesh points in time
    # help variables
    
    Cx2, Cy2 = (dt/dx)**2, (dt/dy)**2
    dt2 = dt**2
    
    # Allow f and V to be None or 0
    if f is None or f == 0:
        f = (lambda x, y, t: 0) if version == 'scalar' else \
            lambda x, y, t: zeros((x.shape[0], y.shape[1]))
        # or simpler: x*y*0
    if V is None or V == 0:
        V = (lambda x, y: 0) if version == 'scalar' else \
            lambda x, y: zeros((x.shape[0], y.shape[1]))


    u   = zeros((Nx+1,Ny+1))   # stores u[n+1]
    u_1 = zeros((Nx+1,Ny+1))   # stores u[n]
    u_2 = zeros((Nx+1,Ny+1))   # stores u[n-1]
    f_a = zeros((Nx+1,Ny+1))   # for compiled loops

    Ix = range(0, u.shape[0])
    Iy = range(0, u.shape[1])
    It = range(0, t.shape[0])

    import time; t0 = time.clock()          # for measuring CPU time

    # Load initial condition into u_1
    if version == 'scalar':
        for i in Ix:
            for j in Iy:
                u_1[i,j] = I(x[i], y[j])
    else: #vectorized version
        u_1[:,:] = I(xv, yv)

    if user_action is not None:
        user_action(u_1, x, xv, y, yv, t, 0)

    # Special formula for first time step
    n = 0
    # First step requires a special formula, use either the scalar 
    # or vectorized version (the impact of more efficient loops than
    # in advance_vectorized is small as this is only one step)
    if version == 'scalar':
        try :
            u = advance_scalar(
                u, u_1, u_2, f, x, y, t, n, Cx2, Cy2, b,
                q, dt2, V, step1=True)
        except Exception as err:
            print(err)

    else:
        try:
            f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
            V_a = V(xv, yv)
            u = advance_vectorized(
                u, u_1, u_2, f, x, y, t, n, Cx2, Cy2, 
                q, dt2, step1=True)
        except Exception as err:
            print(err)

    if user_action is not None:
        user_action(u, x, xv, y, yv, t, 1)

    # Update data structures for next step
    #u_2[:] = u_1;  u_1[:] = u  # safe, but slower
    u_2, u_1, u = u_1, u, u_2

    for n in It[1:-1]:
        if version == 'scalar':
            # use f(x,y,t) function
            u = advance(u, u_1, u_2, f, x, y, t, n, Cx2, Cy2, b,
                q, dt2)
        else:
            f_a[:,:] = f(xv, yv, t[n])  # precompute, size as u
            u = advance(u, u_1, u_2, f_a, Cx2, Cy2, dt2)

        if user_action is not None:
            if user_action(u, x, xv, y, yv, t, n+1):
                break
        # Update data structures for next step
        #u_2[:] = u_1;  u_1[:] = u  # safe, but slower
        u_2, u_1, u = u_1, u, u_2

    # Important to set u = u_1 if u is to be returned!
    t1 = time.clock()
    # dt might be computed in this function so return the value
    #print u
    return dt, t1 - t0
    



def advance_scalar(u, u_1, u_2, f, x, y, t, n, Cx2, Cy2, b,
                   q, dt2, V=True, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])  
    if step1:
        dt = sqrt(dt2)
        Cx2 = 0.25*Cx2*(1+0.5*b*dt); Cy2 = 0.25*Cy2*(1+0.5*b*dt)
        dt2 = 0.5*dt2*(1+0.5*b*dt)
        D1 = 1; D2 = 0
    else:
        dt = sqrt(dt2)
        Cx2 = (1./(1+0.5*b*dt))*Cx2; Cy2 = (1./(1+0.5*b*dt))*Cy2
        dt2 = dt2/(1+0.5*b*dt)
        D1 = 2/(1+0.5*b*dt); D2 = (1-0.5*b*dt)/(1+0.5*b*dt)
    
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            im, jm, ip, jp = i-1, j-1, i+1, j+1
            #im = i+1 if i == 0 else i-1
            #jm = j+1 if j == 0 else j-1
            #ip = i-1 if i == Ix[-1] else i+1
            #jp = j-1 if j == Iy[-1] else j+1
            q_xu_xx = (q(x[i],y[j]) + q(x[ip],y[j]))*(u_1[ip,j]-u_1[i,j]) \
                      -(q(x[i],y[j]) + q(x[im],y[j]))*(u_1[i,j]-u_1[im,j])
            q_yu_yy = (q(x[i],y[j]) + q(x[i],y[jp]))*(u_1[i,jp]-u_1[i,j]) \
                      -(q(x[i],y[j]) + q(x[i],y[jm]))*(u_1[i,j]-u_1[i,jm])
            u[i,j]  = D1*u_1[i,j] - D2*u_2[i,j] + Cx2*q_xu_xx + Cy2*q_yu_yy \
                        + dt2*f(x[i], y[j], t[n])
            if step1:
                u[i,j] += (1-0.5*b*dt)*dt*V(x[i], y[j])
    
    
    # Boundary condition u=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0
    #print("solution", u)
    return u

        
def advance_vectorized(u, u_1, u_2, f, x, y, t, n, Ix, Iy, Cx2, Cy2, q,
                   dt2, A, B, C, C0, C1, C2, V=None, step1=False):
    pass
    

def quadratic(Nx, Ny, version):
    """Exact discrete solution of the scheme."""

    def exact_solution(x, y, t):
        return x*(Lx - x)*y*(Ly - y)*(1 + 0.5*t)

    def I(x, y):
        return exact_solution(x, y, 0)

    def V(x, y):
        return 0.5*exact_solution(x, y, 0)

    def f(x, y, t):
        return 2*c**2*(1 + 0.5*t)*(y*(Ly - y) + x*(Lx - x))

    Lx = 5;  Ly = 2
    c = 1.5
    dt = -1 # use longest possible steps
    T = 18

    def assert_no_error(u, x, xv, y, yv, t, n):
        u_e = exact_solution(xv, yv, t[n])
        diff = abs(u - u_e).max()
        tol = 1E-12
        msg = 'diff=%g, step %d, time=%g' % (diff, n, t[n])
        assert diff < tol, msg

    new_dt, cpu = solver(
        I, V, f, c, Lx, Ly, Nx, Ny, dt, T,
        user_action=assert_no_error, version=version)
    return new_dt, cpu


def test_quadratic():
    # Test a series of meshes where Nx > Ny and Nx < Ny
    versions = 'scalar', 'vectorized', 'cython', 'f77', 'c_cy', 'c_f2py'
    for Nx in range(2, 6, 2):
        for Ny in range(2, 6, 2):
            for version in versions:
                print 'testing', version, 'for %dx%d mesh' % (Nx, Ny)
                quadratic(Nx, Ny, version)


def gaussian(plot_method=2, version='vectorized', save_plot=True):
    """
    Initial Gaussian bell in the middle of the domain.
    plot_method=1 applies mesh function, =2 means surf, =0 means no plot.
    """
    # Clean up plot files
    import glob
    for name in glob.glob('tmp_*.png'):
        os.remove(name)

    Lx = 10
    Ly = 10
    
    def I(x, y):
        """Gaussian peak at (Lx/2, Ly/2)."""
        return exp(-0.5*(x-Lx/2.0)**2 - 0.5*(y-Ly/2.0)**2)
    def q(x,y):
        return 1.
    #if plot_method == 3:
        #from mpl_toolkits.mplot3d import axes3d
        #import matplotlib.pyplot as plt
        #from matplotlib import cm
        #plt.ion()
        #fig = plt.figure()
        #u_surf = None

    def plot_u(u, x, xv, y, yv, t, n):
        if t[n] == 0:
            time.sleep(2)
        if plot_method == 1:
            mesh(x, y, u, title='t=%g' % t[n], zlim=[-1,1],
                 caxis=[-1,1])
        elif plot_method == 2:
            surfc(xv, yv, u, title='t=%g' % t[n], zlim=[-1, 1],
                  colorbar=True, colormap=hot(), caxis=[-1,1],
                  shading='flat')
        elif plot_method == 3:
            fig = plt.figure()
            ax  = fig.add_subplot(111, projection='3d')
            ax.plot_surface(xv, yv, u, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
            #plt.show()

        if plot_method > 0:
            path = 'graphics'
            time.sleep(0) # pause between frames
            if save_plot:
                filename = '{}/tmp_{}04d.png'.format(path,n)
                savefig(filename)  # time consuming!
                

    Nx = 40; Ny = 40; T = 20
    q = lambda x, y: 1
    dt, cpu = solver(I, None, None, q, Lx, Ly, Nx, Ny, -1, T,b=0,
                      user_action=plot_u,version=version)
                      
    #print"time step in gausian", dt
                     
                     


#I, V, f, q, Lx, Ly, Nx, Ny, dt, T, b,
#user_action=None, version='scalar'
if __name__ == '__main__':
    #test_quadratic()
    gaussian(plot_method=2, version='scalar', save_plot=True)
