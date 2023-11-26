
import numpy as np
import sympy as smp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter


xi, k, bw, f0, m, l, t = smp.symbols(r'\xi k \Omega F_0 m l t')


x, f = smp.symbols(r'x f', cls=smp.Function)


# Position
x = x(t)
x_dif = smp.diff(x, t)
x_ddif = smp.diff(x_dif, t)


#Force
f = f(t)
f = f0*smp.sin(bw*t)


w = smp.sqrt(k/m)
c = 2*xi*w*m
eq_struct = smp.Eq(m*x_ddif, -k*(x-l) - c*x_dif + f)




sols = smp.solve(eq_struct, x_ddif, simplify=True, rational=False)



dxdt_f = smp.lambdify(x_dif, x_dif, modules=['numpy']) 
dvdt_f = smp.lambdify((m, k, l, xi, x, x_dif, f0, bw, t), sols[0], modules=['numpy'])


def dSdt(S, t):
    x, v = S
    return [
    dxdt_f(v),
    dvdt_f(m, k, l, xi, x, v, f0, bw, t),
    ]


t = np.linspace(0, 40, 1000)

#Parameters of the structure
m = 8
w = 2
k = np.power(w, 2)*m
l = 10
xi = 0.1

#Force parameters
f0 = 20
bw = 2

#Initial x and v
y0=[10, 0]
ans = odeint(dSdt, y0, t=t)



x = ans.T[0]


#Animation (GIF)
def animate(i):
    ax.clear()
    ax.grid()

    # Define the side length for the square boxes
    box_size = 1.5

    # Draw mass 1
    rect1 = plt.Rectangle((x[i] - box_size / 2, 0), box_size, box_size, edgecolor='black', facecolor='red')
    ax.add_patch(rect1)

    # Draw spring 1
    spring1 = plt.Line2D([0, x[i] - box_size / 2], [box_size / 2, box_size / 2], color='grey', linestyle='-', linewidth=2)
    ax.add_line(spring1)


    # Draw the ground
    ax.add_line(plt.Line2D([0, 14], [0, 0], color='black', linestyle='-', linewidth=2))

    ax.set_ylim(-0.5, 2)
    ax.set_xlim(0, 17)
    ax.set_title(f'Time: {t[i]:.2f}s')

fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ani = animation.FuncAnimation(fig, animate, frames=1000, interval=50)
ani.save('TMD.gif', writer='pillow', fps=25)




