import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from scipy.special import jv, jvp

# Constants 
alpha_degrees = 10
alpha = np.radians(alpha_degrees)
Re_values = [-100, -50, -5, 5, 50, 100]
rho = 1.225
v0 = 1.0
p_infinity = 101325
C_p = 1005
k = 0.0262
nu = v0 / max(np.abs(Re_values))
Pr_number = (rho * C_p * nu) / k

def jeffery_hamel_ode(eta, Y, Re):
    return np.array([Y[1], Y[2], -2 * Re * alpha * Y[0] * Y[1] - 4 * alpha ** 2 * Y[1]])

def boundary_conditions(Ya, Yb):
    return np.array([Ya[0] - 1, Ya[1], Yb[0]])

def pressure_coefficient(eta, Re, alpha):
    F = jv(0, 2 * alpha * eta) / jv(0, 2 * alpha)
    F3 = jvp(0, 2 * alpha * eta, 3) / jv(0, 2 * alpha)
    c_p = 1 + (4 * alpha ** 2) / (alpha * Re) * (1 - F) + 1 / (alpha * Re) * F3
    return c_p

eta = np.linspace(0, 1, 100)

fig_velocity, ax_velocity = plt.subplots(figsize=(8, 4))
fig_pressure, ax_pressure = plt.subplots(figsize=(8, 4))
fig_shear_stress, ax_shear_stress = plt.subplots(figsize=(8, 4))
fig_radial_stress, ax_radial_stress = plt.subplots(figsize=(8, 4))
fig_circumferential_stress, ax_circumferential_stress = plt.subplots(figsize=(8, 4))

for Re in Re_values:
    Y_guess = np.zeros((3, eta.size))
    Y_guess[0] = 1 - eta

    solution = solve_bvp(lambda eta, Y: jeffery_hamel_ode(eta, Y, Re), boundary_conditions, eta, Y_guess)

    if solution.success:
        eta_plot = np.linspace(0, 1, 100)
        Y_plot = solution.sol(eta_plot)

        c_p_values = pressure_coefficient(eta_plot, Re, alpha)
        p_diff = c_p_values * (0.5 * rho * v0 ** 2)
        p_values = p_diff + p_infinity

        tau_theta_r = 2 / Re * Y_plot[1] * (0.5 * rho * v0 ** 2)

        tau_rr = 4 * alpha / Re * Y_plot[0] * (0.5 * rho * v0 ** 2)
        tau_theta_theta = -tau_rr

        ax_velocity.plot(eta_plot, Y_plot[0], label=f'Re={Re}')
        ax_pressure.plot(eta_plot, p_values, label=f'Re={Re}')
        ax_shear_stress.plot(eta_plot, tau_theta_r, label=f'Re={Re}')
        ax_radial_stress.plot(eta_plot, tau_rr, label=f'Re={Re}')
        ax_circumferential_stress.plot(eta_plot, tau_theta_theta, label=f'Re={Re}')
    else:
        print(f'Solution was not successful for Re={Re}. Try adjusting the initial guess or mesh.')

ax_velocity.set_title('Velocity Profile')
ax_velocity.set_xlabel('η')
ax_velocity.set_ylabel('f(η)')
ax_velocity.legend()
ax_velocity.grid(True)

ax_pressure.set_title('Pressure Distribution')
ax_pressure.set_xlabel('η')
ax_pressure.set_ylabel('p')
ax_pressure.legend()
ax_pressure.grid(True)

ax_shear_stress.set_title('Shear Stress')
ax_shear_stress.set_xlabel('η')
ax_shear_stress.set_ylabel('τ_θr')
ax_shear_stress.legend()
ax_shear_stress.grid(True)

ax_radial_stress.set_title('Radial Normal Stress')
ax_radial_stress.set_xlabel('η')
ax_radial_stress.set_ylabel('τ_rr')
ax_radial_stress.legend()
ax_radial_stress.grid(True)

ax_circumferential_stress.set_title('Circumferential Normal Stress')
ax_circumferential_stress.set_xlabel('η')
ax_circumferential_stress.set_ylabel('τ_θθ')
ax_circumferential_stress.legend()
ax_circumferential_stress.grid(True)

fig_velocity.tight_layout()
fig_pressure.tight_layout()
fig_shear_stress.tight_layout()
fig_radial_stress.tight_layout()
fig_circumferential_stress.tight_layout()
plt.show()


def jeffery_hamel_ode_without_dissipation(eta, Y, Re):
    f, df, ddf = Y[:3]
    G, dG = Y[3:5]
    return np.array([df, ddf, -2 * Re * alpha * f * df - 4 * alpha ** 2 * df, dG, (4 + 2 * f * Pr_number) * G])

def jeffery_hamel_ode_with_dissipation(eta, Y, Re):
    f, df, ddf = Y[:3]
    G, dG = Y[3:5]
    return np.array([df, ddf, -2 * Re * alpha * f * df - 4 * alpha ** 2 * df, dG, (4 + 2 * f * Pr_number) * G + Pr_number * (4 * f ** 2 - df ** 2)])

def boundary_conditions(Ya, Yb):
    return np.array([Ya[0] - 1, Ya[1], Yb[0], Ya[3], Yb[3]])

fig_G_without_dissipation, ax_G_without_dissipation = plt.subplots(figsize=(10, 6))

for Re in Re_values:
    Y_guess = np.zeros((5, eta.size))
    Y_guess[0] = 1 - eta

    solution_without_dissipation = solve_bvp(lambda eta, Y: jeffery_hamel_ode_without_dissipation(eta, Y, Re), boundary_conditions, eta, Y_guess)

    if solution_without_dissipation.success:
        eta_plot = np.linspace(0, 1, 100)
        Y_plot = solution_without_dissipation.sol(eta_plot)
        ax_G_without_dissipation.plot(eta_plot, Y_plot[3], label=f'Re={Re}')
    else:
        print(f'Solution without dissipation was not successful for Re={Re}.')

ax_G_without_dissipation.set_xlabel('η')
ax_G_without_dissipation.set_ylabel('G')
ax_G_without_dissipation.set_title('Dimensionless Temperature Function G for Different Reynolds Numbers (Without Dissipation)')
ax_G_without_dissipation.legend()
ax_G_without_dissipation.grid(True)

# Create the plot with dissipation
fig_G_with_dissipation, ax_G_with_dissipation = plt.subplots(figsize=(10, 6))

for Re in Re_values:
    Y_guess = np.zeros((5, eta.size))
    Y_guess[0] = 1 - eta

    solution_with_dissipation = solve_bvp(lambda eta, Y: jeffery_hamel_ode_with_dissipation(eta, Y, Re), boundary_conditions, eta, Y_guess)

    if solution_with_dissipation.success:
        eta_plot = np.linspace(0, 1, 100)
        Y_plot = solution_with_dissipation.sol(eta_plot)
        ax_G_with_dissipation.plot(eta_plot, Y_plot[3], label=f'Re={Re}')
    else:
        print(f'Solution with dissipation was not successful for Re={Re}.')

ax_G_with_dissipation.set_xlabel('η')
ax_G_with_dissipation.set_ylabel('G')
ax_G_with_dissipation.set_title('Dimensionless Temperature Function G for Different Reynolds Numbers (With Dissipation)')
ax_G_with_dissipation.legend()
ax_G_with_dissipation.grid(True)

# Finalize the layout and show the plots
plt.tight_layout()
plt.show()