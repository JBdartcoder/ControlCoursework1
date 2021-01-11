import sympy as sym
import control as C
import matplotlib.pyplot as plt
import numpy as np


# Define all involved symbolic variables and constants
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
def evaluate_constants(f):
    """
    This function subs in
    the constant values into
    A1, A2, A3, B1, B2 and B3
    """
    return f.subs([(m, m_value), (R, R_value), (L_0, L_0_value), (L_1, L_1_value), (alpha, alpha_value),
                   (c, c_value), (delta, delta_value), (k, k_value), (b, b_value), (x3_eq, x3_eq_value),
                   (x1_eq, x1_eq_value), (V_e, V_e_value)])


A_1_value = float(evaluate_constants(A_1))
A_2_value = float(evaluate_constants(A_2))
A_3_value = float(evaluate_constants(A_3))
B_1_value = float(evaluate_constants(B_1))
B_2_value = float(evaluate_constants(B_2))
B_3_value = float(evaluate_constants(B_3))


# Use A_1 -> B_3 to determine the coefficients of the numerator
# and the denominator (from s^3 - s^0 term)
num_value = A_3_value * B_1_value
s_3_den_value = 1.0       # coeff. of s^3
s_2_den_value = -B_3_value - A_2_value      # coeff. of s^2
s_1_den_value = (B_3_value * A_2_value) - A_1_value     # coeff. of s^1
s_0_den_value = (B_3_value*A_1_value) - (A_3_value*B_2_value)       # coeff. of s^0

# Declare overall numerator and denominator of transfer function
num_Gx = [num_value]
den_Gx = [s_3_den_value, s_2_den_value, s_1_den_value, s_0_den_value]

# Declare additional symbols from the transfer function
s, t = sym.symbols('s, t')

# Declare the system as a transfer function signal using the numerator and denominator from above
G_1 = C.TransferFunction(num_Gx, den_Gx)

# Declare the transfer function for the sensor - first order system
# G_2 = kappa / ((tau * s) + 1)
G_2 = C.TransferFunction([1], [tau, 1])


def pid(kp, ki, kd):
    # This function constructs the transfer function of a PID
    # controller with given parameters
    diff = C.TransferFunction([1, 0], 1)
    intgr = C.TransferFunction(1, [1, 0])
    pid_tf = kp + kd * diff + ki * intgr
    return pid_tf


# Setting the PID values
Kp = 70
Ki = 0.1
Kd1 = 3
Kd2 = 4
Kd3 = 5
controller1 = pid(Kp, Ki, Kd1)
controller2 = pid(Kp, Ki, Kd2)
controller3 = pid(Kp, Ki, Kd3)

# Declaring the simulation range and num_points
t_final = 1
num_points = 500
t_span = np.linspace(0, t_final, num_points)

# Using the feedback function to calculate the transfer function of the entire system
G_x = C.feedback(G_1, G_2)

# The transfer functions for the different values of kd
G_d1 = C.feedback(G_x, controller1)
G_d2 = C.feedback(G_x, controller2)
G_d3 = C.feedback(G_x, controller3)

# Using the impulse response function to simulation the controller
t_imp, x_imp1 = C.impulse_response(G_d1, t_span)
_, x_imp2 = C.impulse_response(G_d2, t_span)
_, x_imp3 = C.impulse_response(G_d3, t_span)

# Plotting the graphs
plt.plot(t_imp, x_imp1, label='$K_d = 3$')
plt.plot(t_imp, x_imp2, label='$K_d = 4$')
plt.plot(t_imp, x_imp3, label='$K_d = 5$')
plt.xlabel('Time (s)')
plt.ylabel('$x_1$ (m)')
plt.legend()
plt.grid()
plt.show()
