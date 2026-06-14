import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hermite

#Parameter at beam waist
w_0 = 50.0  # beam waist radius
A = 1.0

#Inputs
m, n = map(int, input("Enter mode indices (m n): ").split())

#Grid
x = np.linspace(-4*w_0, 4*w_0, 500)
y = np.linspace(-4*w_0, 4*w_0, 500)
X,Y = np.meshgrid(x,y)


#Hermite Polynomials
hm = hermite(m)
hn = hermite(n)

#Electric Field 
E = A * hm(np.sqrt(2)*X/w_0) * hn(np.sqrt(2)*Y/w_0) *  np.exp(-(X**2 + Y**2)/w_0**2)

#Intensity
I = np.abs(E)**2

#plot
plt.imshow(I, cmap='inferno',
           extent=[x.min(), x.max(), y.min(), y.max()],
           origin='lower')
plt.colorbar(label='Intensity')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'TEM({m},{n})')
plt.show()
