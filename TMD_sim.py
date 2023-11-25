
import numpy as np
import sympy as smp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter

# Symbolic Parameters
xi, xi_d, k, k_d, l1, l2, bw, f0, m, m_d, t = smp.symbols(r'\xi \xi_d k k_d l_1 l_2 \Omega F_0 m m_d t')
u, ud, x1, x2, f = smp.symbols(r'u u_d x_1 x_2 f', cls=smp.Function)


# U
u = u(t)
u_dif = smp.diff(u, t)
u_ddif = smp.diff(u_dif, t)


#U_d
ud = ud(t)
ud_dif = smp.diff(ud, t)
ud_ddif = smp.diff(ud_dif, t)

#Positions

x1 = x1(u)
x2 = x2(u, ud)
x1 = u
x2 = u + ud

#Force
f = f(t)
f = f0*smp.sin(bw*t)


mbar = m_d/m
w = smp.sqrt(k/m)
w_d = smp.sqrt(k_d/m_d)
c = 2*xi*w*m
c_d = 2*xi_d*w_d*m_d


# Movement equations
eq_struct = smp.Eq(m*u_ddif,  -k*(u - l1) - c*u_dif + k_d*(ud - l2) + c_d*ud_dif + f)
eq_tunned = smp.Eq(m_d*(u_ddif + ud_ddif),  - k_d*(ud - l2) - c_d*ud_dif)



sols = smp.solve([eq_struct, eq_tunned], (u_ddif, ud_ddif), simplify = True, rational= False)



dudt_f = smp.lambdify(u_dif, u_dif, modules=['numpy']) 
dvdt_f = smp.lambdify((xi, xi_d, k, k_d, l1, l2, bw, f0, m, m_d, u, u_dif, ud, ud_dif, t), sols[u_ddif], modules=['numpy']) 

du_ddt_f = smp.lambdify(ud_dif, ud_dif, modules=['numpy'])
dv_ddt_f = smp.lambdify((xi, xi_d, k, k_d, l1, l2, bw, f0, m, m_d, u, u_dif, ud, ud_dif, t), sols[ud_ddif], modules=['numpy'])


def dSdt(S, t):
    u, v, ud, vd = S
    return [
    dudt_f(v),
    dvdt_f(xi, xi_d, k, k_d, l1, l2, bw, f0, m, m_d, u, v, ud, vd, t),
    du_ddt_f(vd),
    dv_ddt_f(xi, xi_d, k, k_d, l1, l2, bw, f0, m, m_d, u, v, ud, vd, t)
    ]


t = np.linspace(0, 40, 1000)

#Parameters of the structure
m = 8
w = 2
k = np.power(w, 2)*m
l1 = 10
xi = 0.1

#Parameters of tunned mass
m_d= 4
w_d = 2
k_d = np.power(w_d, 2)*m_d
l2 = 10
xi_d = 0.2

#Force parameters
f0 = 20
bw = 2

#Intial u, v, u_d and v_d
y0=[10, 0, 10 ,0]
ans = odeint(dSdt, y0, t=t)



u = ans.T[0]
ud = ans.T[2]


def get_x1x2(u, ud):
    return(u, u+ud)

x1, x2 = get_x1x2(u, ud)


#Amplitudes
A1 = x1-y0[0]
A2 = x2-y0[0]-y0[2]


#Animation (GIF)
def animate(i):
    ax.clear()
    ax.grid()

    # Define the side length for the square boxes
    box_size = 1.5

    # Draw mass 1
    rect1 = plt.Rectangle((x1[i] - box_size / 2, 0), box_size, box_size, edgecolor='black', facecolor='red')
    ax.add_patch(rect1)

    # Draw spring 1
    spring1 = plt.Line2D([0, x1[i] - box_size / 2], [box_size / 2, box_size / 2], color='grey', linestyle='-', linewidth=2)
    ax.add_line(spring1)

    # Draw mass 2
    rect2 = plt.Rectangle((x2[i] - box_size / 2, 0), box_size, box_size, edgecolor='blue', facecolor='blue')
    ax.add_patch(rect2)

    # Draw spring 2
    spring2 = plt.Line2D([x1[i] + box_size / 2, x2[i] - box_size / 2], [box_size / 2, box_size / 2], color='grey', linestyle='-', linewidth=2)
    ax.add_line(spring2)

    # Draw the ground
    ax.add_line(plt.Line2D([0, 22], [0, 0], color='black', linestyle='-', linewidth=2))

    ax.set_ylim(-0.5, 2)
    ax.set_xlim(0, 25)
    ax.set_title(f'Time: {t[i]:.2f}s')

fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ani = animation.FuncAnimation(fig, animate, frames=1000, interval=50)
ani.save('TMD.gif', writer='pillow', fps=25)





