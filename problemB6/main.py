import sympy as sym
import control as C
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


"""
# Define all involved symbolic variables
# constants
m, g, d, delta, r, R, L_0, L_1, alpha, c, k, b, phi, tau, kappa = sym.symbols('m, g, d, delta, r, R, L_0, L_1, alpha, c, k, b, phi, tau, kappa')
# system variables
x1_eq, x2_eq, x3_eq, V_e = sym.symbols('x1_eq, x2_eq, x3_eq, V_e')

# Declare determined values of A1 -> B3
A_1 = (5 / (7 * m)) * ((2 * c * (x3_eq ** 2) / (delta - x1_eq) ** 3) - k)
A_2 = -(5 * b) / (7 * m)
A_3 = (10 * c * x3_eq) / (7 * m * ((delta - x1_eq) ** 2))
B_1 = 1 / (L_0 + (L_1 * sym.exp(-alpha * (delta - x1_eq))))
B_2 = -(L_1 * alpha * ((-R * x3_eq) + V_e) * sym.exp(-alpha * (delta - x1_eq))) / (
        (L_0 + (L_1 * sym.exp(-alpha * (delta - x1_eq)))) ** 2)
B_3 = -R / (L_0 + (L_1 * sym.exp(-alpha * (delta - x1_eq))))

# State values of given parameters
m_value = 0.425
g_value = 9.81
d_value = 0.42
delta_value = 0.65
r_value = 0.125
R_value = 53
L_0_value = 0.12
L_1_value = 0.25
alpha_value = 1.2
c_value = 6.815
k_value = 1880
b_value = 10.4
phi_value = 42
V_e_value = 36.04
tau = 0.03

# State equilibrium values as determined in derivations
x1_eq_value = 0.47861
x2_eq_value = 0
x3_eq_eqn = V_e / R
x3_eq_value = x3_eq_eqn.subs([(V_e, V_e_value), (R, R_value)])

# Substitute values of the constants into the equations for A_1 -> B_3
A_1_value = A_1.subs(
    [(m, m_value), (c, c_value), (x3_eq, x3_eq_value), (delta, delta_value), (x1_eq, x1_eq_value), (k, k_value)])
A_2_value = A_2.subs([(b, b_value), (m, m_value)])
A_3_value = A_3.subs([(c, c_value), (m, m_value), (x3_eq, x3_eq_value), (delta, delta_value), (x1_eq, x1_eq_value)])
B_1_value = B_1.subs(
    [(L_0, L_0_value), (L_1, L_1_value), (alpha, alpha_value), (delta, delta_value), (x1_eq, x1_eq_value)])
B_2_value = B_2.subs(
    [(L_0, L_0_value), (L_1, L_1_value), (alpha, alpha_value), (delta, delta_value), (x1_eq, x1_eq_value), (R, R_value),
     (x3_eq, x3_eq_value), (V_e, V_e_value)])
B_3_value = B_3.subs(
    [(L_0, L_0_value), (L_1, L_1_value), (alpha, alpha_value), (delta, delta_value), (x1_eq, x1_eq_value),
     (R, R_value)])

# Declare additional symbols from the transfer function
s, t = sym.symbols('s, t')
w = sym.symbols('w', real=True)  # w represents omega
# a-d must be positive for inverse Laplace transform to function correctly
a, b, c, d = sym.symbols('a:d', real=True, positive=True)

# Declare transfer functions G_theta and G_x
# from derivation on additional notes 2
G_1 = (A_3_value * B_1_value) / (((s**2 - (s*A_2_value) - A_1_value) * (s - B_3_value)) - (A_3_value * B_2_value))
# sym.pprint(G_1)

G_2 = kappa / ((tau * s) + 1)
# sym.pprint(G_2)

G_x = G_1 / (1 + (G_1 * G_2))
sym.pprint(G_x.simplify())
# print(sym.latex(G_x.simplify()))
"""


# Define all involved symbolic variables
# constants
m, g, d, delta, r, R, L_0, L_1, alpha, c, k, b, phi, tau, kappa = sym.symbols('m, g, d, delta, r, R, L_0, L_1,'
                                                                              ' alpha, c, k, b, phi, tau, kappa')
# system variables
x1_eq, x2_eq, x3_eq, V_e = sym.symbols('x1_eq, x2_eq, x3_eq, V_e')

# Declare determined values of A1 -> B3
A_1 = (5 / (7 * m)) * ((2 * c * (x3_eq ** 2) / (delta - x1_eq) ** 3) - k)
A_2 = -(5 * b) / (7 * m)
A_3 = (10 * c * x3_eq) / (7 * m * ((delta - x1_eq) ** 2))
B_1 = 1 / (L_0 + (L_1 * sym.exp(-alpha * (delta - x1_eq))))
B_2 = -(L_1 * alpha * ((-R * x3_eq) + V_e) * sym.exp(-alpha * (delta - x1_eq))) / (
        (L_0 + (L_1 * sym.exp(-alpha * (delta - x1_eq)))) ** 2)
B_3 = -R / (L_0 + (L_1 * sym.exp(-alpha * (delta - x1_eq))))

# State values of given parameters
m_value = 0.425
g_value = 9.81
d_value = 0.42
delta_value = 0.65
r_value = 0.125
R_value = 53
L_0_value = 0.12
L_1_value = 0.25
alpha_value = 1.2
c_value = 6.815
k_value = 1880
b_value = 10.4
phi_value = 42
V_e_value = 36.04
tau = 0.03

# State equilibrium values as determined in derivations
x1_eq_value = 0.47861
x2_eq_value = 0
x3_eq_eqn = V_e / R
x3_eq_value = x3_eq_eqn.subs([(V_e, V_e_value), (R, R_value)])

# Substitute values of the constants into the equations for A_1 -> B_3
A_1_value = A_1.subs(
    [(m, m_value), (c, c_value), (x3_eq, x3_eq_value), (delta, delta_value), (x1_eq, x1_eq_value), (k, k_value)])
A_2_value = A_2.subs([(b, b_value), (m, m_value)])
A_3_value = A_3.subs([(c, c_value), (m, m_value), (x3_eq, x3_eq_value), (delta, delta_value), (x1_eq, x1_eq_value)])
B_1_value = B_1.subs(
    [(L_0, L_0_value), (L_1, L_1_value), (alpha, alpha_value), (delta, delta_value), (x1_eq, x1_eq_value)])
B_2_value = B_2.subs(
    [(L_0, L_0_value), (L_1, L_1_value), (alpha, alpha_value), (delta, delta_value), (x1_eq, x1_eq_value), (R, R_value),
     (x3_eq, x3_eq_value), (V_e, V_e_value)])
B_3_value = B_3.subs(
    [(L_0, L_0_value), (L_1, L_1_value), (alpha, alpha_value), (delta, delta_value), (x1_eq, x1_eq_value),
     (R, R_value)])

# Use A_1 -> B_3 to determine the coefficients of the numerator
# and the denominator (from s^3 - s^0 term)

num_value = A_3_value * B_1_value
s_3_den_value = 1.0       # coeff. of s^3
s_2_den_value = -B_3_value - A_2_value      # coeff. of s^2
s_1_den_value = (B_3_value * A_2_value) - A_1_value     # coeff. of s^1
s_0_den_value = (B_3_value*A_1_value) - (A_3_value*B_2_value)       # coeff. of s^0

"""
print(num_value)
print(s_3_den_value)
print(s_2_den_value)
print(s_1_den_value)
print(s_0_den_value)
"""

"""
# Declare overall numerator and denominator of transfer function
num_Gx = [num_value]
den_Gx = [s_3_den_value, s_2_den_value, s_1_den_value, s_0_den_value]
"""
num_Gx = [1639.10919772563]
den_Gx = [1., 181.298987845931, 3919.09863700501, 172943.180855727]

# Declare additional symbols from the transfer function
s, t = sym.symbols('s, t')

# Declare the system as a transfer function signal using the numerator and denominator from above
G_1 = C.TransferFunction(num_Gx, den_Gx)

# Declare the transfer function for the sensor - first order system
# G_2 = kappa / ((tau * s) + 1)
G_2 = C.TransferFunction([1], [tau, 1])

# Combining the two transfer function in parallel
G_x = C.parallel(G_1, G_2)
# G_x = G_1 / (1 + (G_1 * G_2))


def pid(kp, ki, kd):
    # This function constructs the transfer function of a PID
    # controller with given parameters
    diff = C.TransferFunction([1, 0], 1)
    intgr = C.TransferFunction(1, [1, 0])
    pid_tf = kp + kd * diff + ki * intgr
    return pid_tf


Kp = 0.01
Ki = 0.01
Kd = 0.01
controller = -pid(Kp, Ki, Kd)

t_final = 1
num_points = 500
t_span = np.linspace(0, t_final, num_points)


G_d = C.feedback(G_x, controller)

t_imp, x_imp = C.impulse_response(G_d, t_span)


plt.plot(t_imp, x_imp)
plt.xlabel('Time (s)')
plt.ylabel('')
plt.grid()
plt.show()
