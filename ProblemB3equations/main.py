import sympy as sym

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
G_x = (A_3_value * B_1_value) / (((s**2 - (s*A_2_value) - A_1_value) * (s - B_3_value)) - (A_3_value * B_2_value))

# Perform an impulse (kick), step (push)
# Kick -> Dirac pulse: F_s = 1
# Push -> Step response: F_s = 1/s
F_s_impulse = 1
F_s_step = 1 / s

# G_theta responses for kick, push
x_s_impulse = G_x * F_s_impulse
x_t_impulse = sym.inverse_laplace_transform(x_s_impulse, s, t)

x_s_step = G_x * F_s_step
x_t_step = sym.inverse_laplace_transform(x_s_step, s, t)

# Print impulse/step response equations inc. LaTex code
print('x response - impulse')
sym.pprint(x_t_impulse.simplify())
print(sym.latex(x_t_impulse.simplify()))
print('x response - step')
sym.pprint(x_t_step.simplify())
print(sym.latex(x_t_step.simplify()))