import sympy as sym
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

# Define all involved symbolic variables
# constants
m, g, d, delta, r, R, L_0, L_1, alpha, c, k, b, phi = sym.symbols('m, g, d, delta, r, R, L_0, L_1, alpha, c, k, b, phi')
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

# State equilibrium values as determined in derivations
x1_eq_value = 0.47861
x2_eq_value = 0
x3_eq_eqn = V_e / R
x3_eq_value = x3_eq_eqn.subs([(V_e, V_e_value), (R, R_value)])
#x3_eq_value = 0.68

# Substitute values of the constants into the equations for A_1 -> B_3
A_1_value = float(A_1.subs(
    [(m, m_value), (c, c_value), (x3_eq, x3_eq_value), (delta, delta_value), (x1_eq, x1_eq_value), (k, k_value)]))
A_2_value = float(A_2.subs([(b, b_value), (m, m_value)]))
A_3_value = float(A_3.subs([(c, c_value), (m, m_value), (x3_eq, x3_eq_value), (delta, delta_value), (x1_eq, x1_eq_value)]))
B_1_value = float(B_1.subs(
    [(L_0, L_0_value), (L_1, L_1_value), (alpha, alpha_value), (delta, delta_value), (x1_eq, x1_eq_value)]))
B_2_value = float(B_2.subs(
    [(L_0, L_0_value), (L_1, L_1_value), (alpha, alpha_value), (delta, delta_value), (x1_eq, x1_eq_value), (R, R_value),
     (x3_eq, x3_eq_value), (V_e, V_e_value)]))
B_3_value = float(B_3.subs(
    [(L_0, L_0_value), (L_1, L_1_value), (alpha, alpha_value), (delta, delta_value), (x1_eq, x1_eq_value),
     (R, R_value)]))
print(B_3_value)

# Use A_1 -> B_3 to determine the coefficients of the numerator
# and the denominator (from s^3 - s^0 term)

num_value = A_3_value * B_1_value
s_3_den_value = 1.       # coeff. of s^3
s_2_den_value = -B_3_value - A_2_value      # coeff. of s^2
s_1_den_value = (B_3_value * A_2_value) - A_1_value     # coeff. of s^1
s_0_den_value = (B_3_value*A_1_value) - (A_3_value*B_2_value)       # coeff. of s^0
num_Gx = [num_value]
den_Gx = [s_3_den_value, s_2_den_value, s_1_den_value, s_0_den_value]
G_theta_tf = signal.TransferFunction(num_Gx, den_Gx)
t, y = signal.step(G_theta_tf)
plt.figure(1)
plt.plot(t,y,'r-')
plt.show()

t,y = signal.impulse(G_theta_tf)
plt.figure(2)
plt.plot(t,y,'r-')
plt.show()
# F_s_kick = 1
# F_s_push = 1/s
# # Determine response of G_theta
# t_step, y_step, x_step = ctrl.forced_response(G_theta_tf, t_span, F_s_push)



#
# # Declare overall numerator and denominator of transfer function
# num_Gx = [num_value]
# den_Gx = [s_3_den_value, s_2_den_value, s_1_den_value, s_0_den_value]
#
# # Declare the system as a transfer function signal using the numerator and denominator from above
# sys = signal.TransferFunction(num_Gx, den_Gx)

# # Plot the transfer function impulse response against time
# t = np.linspace(0, 1, 1000)   # Up to 1s, 1000 intervals
# t, y = signal.impulse(sys, T=t)
# plt.figure(1)
# plt.plot(t,y,'k-',label="impulse response")
# plt.xlabel('t (seconds)')
# plt.ylabel('y(t) (unit unknown)')
# plt.legend(loc='best')
# plt.show()
#
# # Plot the transfer function step response against time
# t, y = signal.step(sys, T=t)
# plt.figure(1)
# plt.plot(t,y,'k-',label="step response")
# plt.xlabel('t (seconds)')
# plt.ylabel('y(t) (unit unknown)')
# plt.legend(loc='best')
# plt.show()
#
