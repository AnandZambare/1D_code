import matplotlib.pyplot as plt
from numpy import *
from celluloid import Camera
# Computer assignment 1 Foundations of CFD
# Roll no - AM21S004
# Explicit FTCS method for unsteady heat conduction equation with Dirichlet boundary conditions

# Domain
L = float(input('Enter the length of domain'))                            # Length of domain from user
Nx = int(input('Enter the no. of grid points in the domain'))           # No. of grid points
alpha = 1                                                          # Thermal diffusivity
dt = float(input('Enter the time step size'))                             # time step size from user
dx = L/(Nx-1)                                                      # delta X in the domain
# using delta X form a grid
X = arange(0, L+dx, dx)
# Temperatures for corresponding grid points
T = zeros(Nx)                                   # This is also the given initial condition
gamma = (alpha*dt/(dx**2))                     # constant to be used in equation
time = 0                                       # initialize time
error = ones(Nx)                               # Error vector initialize
fig = plt.figure()
camera = Camera(fig)
plt.plot(X, T, color='RED')
err = average(error)


# Time Marching
while err > (10**-6):
    time = time + dt
    T[0] = 1.0  # Type 1 boundary condition
    T_prev = T.copy()
    for i in range(1, Nx-1):
        T[i] = (gamma * T_prev[i + 1]) + ((1 - (2 * gamma)) * T_prev[i]) + (gamma * T_prev[i - 1])
    T[0] = 1.0  # Type 1 boundary condition
    T[Nx - 1] = 0.0  # Type 1 boundary condition
    for i in range(Nx):
        error[i] = T[i] - T_prev[i]
    plt.plot(X, T, color='RED')
    camera.snap()
    err = average(error)
animation = camera.animate(interval=1)
fig = plt.figure()
# plt.show()
plt.grid()
plt.plot(X, T)
# plt.show()
print( 'Steady state temperature values are', T)