import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

alpha_degrees = 10
alpha = np.radians(alpha_degrees)
Re_values = [-100, -50, -5, 5, 50, 100]
rho = 1.225
v0 = 1.0
C_p = 1005
k = 0.0262
nu = v0 / max(np.abs(Re_values))
Pr_number = (rho * C_p * nu) / k
Ts = 300

def jeffery_hamel_ode(eta, Y, Re):
    f, df, ddf = Y[:3]
    G, dG = Y[3:5]
    return np.array([df, ddf, -2 * Re * alpha * f * df - 4 * alpha ** 2 * df, dG, (4 + 2 * f * Pr_number) * G + Pr_number * (4 * f ** 2 - df ** 2)])

def boundary_conditions(Ya, Yb):
    return np.array([Ya[0] - 1, Ya[1], Yb[0], Ya[3], Yb[3]])

eta = np.linspace(0, 1, 100)

fig_G, ax_G = plt.subplots(figsize=(8, 4))

for Re in Re_values:
    Y_guess = np.zeros((5, eta.size))
    Y_guess[0] = 1 - eta

    solution = solve_bvp(lambda eta, Y: jeffery_hamel_ode(eta, Y, Re), boundary_conditions, eta, Y_guess)

    if solution.success:
        eta_plot = np.linspace(0, 1, 100)
        Y_plot = solution.sol(eta_plot)

        ax_G.plot(eta_plot, Y_plot[3], label=f'Re = {Re}')

ax_G.set_xlabel('η')
ax_G.set_ylabel('G')
ax_G.set_title('Dimensionless Temperature Function G for Different Reynolds Numbers')
ax_G.legend()
ax_G.grid(True)

plt.show()
