import sympy as sym
import matplotlib.pyplot as plt
from scipy import signal as sig

# Define all involved symbolic variables
# constants
m, g, d, delta, r, R, L_0, L_1, alpha, c, k, b, phi = sym.symbols('m, g, d, delta, r, R, L_0, L_1, alpha, c, k, b, phi')
# system variables
x1_eq, x2_eq, x3_eq, V_e = sym.symbols('x1_eq, x2_eq, x3_eq, V_e')

# Declare determined values of A1 -> B3
# As determined in Part A
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
V_e_value = 36.0588

# State equilibrium values as determined in derivations
x1_eq_value = 0.47861  # Previously calculated in Part B
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
s_3_den_value = 1.  # coeff. of s^3
s_2_den_value = -B_3_value - A_2_value  # coeff. of s^2
s_1_den_value = (B_3_value * A_2_value) - A_1_value  # coeff. of s^1
s_0_den_value = (B_3_value * A_1_value) - (A_3_value * B_2_value)  # coeff. of s^0
num_Gx = [num_value]
den_Gx = [s_3_den_value, s_2_den_value, s_1_den_value, s_0_den_value]

# Declare overall transfer function
G_x_tf = sig.TransferFunction(num_Gx, den_Gx)

# # Use scipy signal library to plot impulse response
t, y = sig.impulse(G_x_tf)
plt.figure(2)
plt.plot(t, y, 'r-')
plt.title('Impulse response of transfer function')
plt.xlabel('time (s)')
plt.ylabel('$x_{1}$ (m)')
plt.grid()
plt.show()

# Use scipy signal library to plot step response
t, y = sig.step(G_x_tf)
plt.figure(1)
plt.plot(t, y, 'g-')
plt.title('Step response of transfer function')
plt.xlabel('time (s)')
plt.ylabel('$x_{1}$ (m)')
plt.grid()
plt.show()
