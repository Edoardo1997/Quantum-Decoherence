import numpy as np
import  matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import scipy.sparse
from scipy.sparse.linalg import inv
import subprocess
from termcolor import colored
import warnings
import pickle
from numba import jit   
import os

# Settings

omega = 1.  # Harmonic Oscillator frequency

precision = 50  
# Precision should be multiple of 50, I used 150 
# in serious simulations, 50 can be used for a fast run

gamma = 0.001  # Underdamped regime
beta = 0.01  # kT >> hbar omega

dt = 0.0001
dx = 0.01
alpha1 = -0.25j*dt/dx**2
beta = 0.25*gamma*dt/dx
N = 100  # In serious simulation I used 500, 100 can be used for fast runs
num = 20 


# Function used to initialize the system


def mu(s):
    """
    Return \mu from squeezing parameter
    """
    beta=(np.exp(2*s)+np.exp(-2*s))/2.
    return (beta-1.)/(beta+1)


def alpha(s):
    """
    Return (1-mu)/(1+mu)
    """
    return 2./(np.exp(2*s)+np.exp(-2*s))


def squeezed(squeezing, x):
    """
    Return normalized squeezed state of squeezing s
    """
    a=alpha(squeezing)
    return np.exp(-(x**2*omega/2*a))*(a*omega/np.pi)**(1./4)


# Function to process the data produced

def integrate(x, y):
    """
    Basic integration of y(x) over x
    """
    dx = x[1:]-x[:-1]
    y_m = (y[1:]+y[:-1])/2
    return np.dot(dx, y_m)


def trace(x,matrix):
    """
    Trace for continuous operators O(x)
    """
    return integrate(x, matrix.diagonal())


def product(x, matrix1, matrix2):
    """
    Product for continuous matrices A(x) B(x)
    """
    dx = x[1:]-x[:-1]
    m1_m = (matrix1[:, 1:]+matrix1[:, :-1])/2
    m2_m = (matrix2[1:, :]+matrix2[:-1, :])/2
    return np.einsum('ij,jk,j->ik', m1_m, m2_m, dx)


def tr2(x, squeezings):
    """
    Used for debugging purpuses, check that the inizial states are normalized.
    """
    L = []
    for s in squeezings:
        Dx = np.exp(s)/2
        lim = 10*Dx
        x = np.linspace(-lim, lim, precision)
        s_state = squeezed(s, x)
        s_matrix = np.outer(s_state, s_state)
        L.append(trace(x, product(x, s_matrix, s_matrix)))
    L = np.asarray(L)
    return L


# Functions for numerical evolution of the system: [algorithm Crank-Nicolson ]


def V_tilde(x, y):
    """
    Return the potential acting on the system
    """
    X, Y = np.meshgrid(x, y)
    X, Y = X.T, Y.T
    return 1./2*omega*(X**2-Y**2)-1j*gamma/beta*(X-Y)**2


def exp_V(dt, x, y):
    """
    Return the propagator of the system
    """
    return np.exp(-0.5j*dt*V_tilde(x, y))


def D2_x(rho, dx):
    """
    Second derivative on x for Crank-Nicolson
    """
    return (np.roll(rho, -1, axis=0)+np.roll(rho, 1, axis=0)-2*rho)/dx**2


def D2_y(rho, dx):
    """
    Second derivative on y for Crank-Nicolson
    """
    return (np.roll(rho, -1, axis=1)+np.roll(rho, 1, axis=1)-2*rho)/dx**2


def D_x(rho, dx):
    """
    First derivative on x for Crank-Nicolson
    """
    return (np.roll(rho, -1, axis=0)-np.roll(rho, 1, axis=0))/(2*dx)


def D_y(rho, dx):
    """
    First derivative on y for Crank-Nicolson
    """
    return (np.roll(rho, -1, axis=1)-np.roll(rho, 1, axis=1))/(2*dx)


def L_x(rho, X, Y, dx):
    """
    Derivative on x of the Limdladian operator
    """
    return 0.5j*D2_x(rho, dx)-0.5*gamma*(X-Y)*D_x(rho, dx)


def L_y(rho, X, Y, dx):
    """
    Derivative on y of the Limdladian operator
    """
    return -0.5j*D2_y(rho, dx)+0.5*gamma*(X-Y)*D_y(rho, dx)


@jit
def Thomas(l,m,u,sol):
    """
    Implementation of thomas algorithm to solve tridiagonal sistems
    l = Lower diagonal,
    m = Main Diag,
    u = Upper Diag,
    sol = solution vector
    Tested on random matrix, work fine, much more effective than Gauss
    Complessity O(N)
    """
    n = len(sol)
    w = np.zeros(n-1, float)
    g = np.zeros(n, float)
    p = np.zeros(n,float)

    w[0] = u[0]/m[0]
    g[0] = sol[0]/m[0]

    for i in range(1,n-1):
        w[i] = u[i]/(m[i] - l[i-1]*w[i-1])
    for i in range(1,n):
        g[i] = (sol[i] - l[i-1]*g[i-1])/(m[i] - l[i-1]*w[i-1])
    p[n-1] = g[n-1]
    for i in range(n-1,0,-1):
        p[i-1] = g[i-1] - w[i-1]*p[i]
    return p


def diagonals(alpha, beta, x, y):
    """
    Prepare the diagonal that must be passed to Thomson algorithm
    """
    down = alpha-beta*(x-y)
    dmain = (1-2*alpha)*np.ones(len(x))
    dup = alpha+beta*(x-y)
    return down, dmain, dup


def timestep(rho, X, Y, dx, dt):
    """
    One single step of the Crank-Nicolson algorithm
    """
    expV = exp_V(dt, X[:, 0], Y[0, :])
    xs = X[:, 0]
    K = rho+0.5*dt*L_y(rho, X, Y, dx)
    alpha1 = -0.25j*dt/dx**2
    beta = 0.25*gamma*dt/dx
    L = []
    for y in range(len(Y[0])):
        down, dmain, dup = diagonals(alpha1, beta, X[:, y], Y[:, y])
        L.append(Thomas(down, dmain, dup, K[:, y]))
    K0 = np.asarray(L).T
    K01 = expV*K0
    K1 = K01+0.5*dt*L_x(K01, X, Y, dx)
    down, dmain, dup = diagonals(-alpha1,-beta,X,Y)  # The change of sign in alpha and beta follows from analitical calculation
    L = []
    for x in range(len(X[:, 0])):
        down, dmain, dup = diagonals(-alpha1, -beta, X[x, :], Y[x, :])  # The change of sign in alpha and beta follows from analitical calculation
        L.append(Thomas(down, dmain, dup, K1[x, :]))
    K00 = np.asarray(L)
    K2 = expV*K00
    T1 = trace(xs, K2)
    T2 = trace(xs,product(xs,K2,K2))
    return K2, T1, T2


def cycle(rho, X, Y, dx, dt, N):
    """
    Complete Crank-Nicolson algorithm
    """
    # Initialization
    rho_next, t1, t2 = timestep(rho, X, Y, dx, dt)
    T1 = [t1]
    T2 = [t2]
    # Iteration of single step
    for i in range(N):
        rho_next, t1, t2 = timestep(rho_next, X, Y, dx, dt)
        T1.append(t1)
        T2.append(t2)
        # Visualizing the status of the computation (since it can take some time)
        if(i%(N/10)==0 and i!=0):
            n = int(i/N*10)
            print('['+n*colored('*', 'red')+(10-n)*'*'+']')
    print('['+10*colored('*', 'red')+']')
    return np.asarray(T1), np.asarray(T2)


x_db = np.array([i for i in range(N+1)])*dt  # used for debugging, now useless

warnings.filterwarnings("ignore") #ignore warnings (dynamical visualization otherwise can be annoying)

solutions = []
squeezings = np.linspace(0, 5, num)  # Inital states

# Iterate entire Crank-Nicolson algorithm over all the initial states
for i, s in enumerate(squeezings):
    Dx = np.exp(s)/2
    limit = 10*Dx
    x = np.linspace(-limit, limit, precision)
    y = np.linspace(-limit, limit, precision)
    s_state = squeezed(s, x)
    s_matrix = np.outer(s_state, s_state)
    X, Y = np.meshgrid(x, y)
    X, Y = X.T, Y.T
    dx = x[1]-x[0]
    T1, T2 = cycle(s_matrix, X, Y, dx, dt, N)
    solutions.append([T1, T2]) # Evaluate trace and trace(rho**2), the second one indicate effect of decoherence
    print(str(i+1)+'/'+str(num))

# Saving data on pikle output file

base = os.getcwd()
path = base + '/pickle_out.txt'

out_file = open(path, 'wb+')
pickle.dump([x_db, solutions], out_file)
out_file.close()

