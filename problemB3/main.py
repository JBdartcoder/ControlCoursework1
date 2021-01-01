import sympy as sym
import numpy as np

# Define all involved symbolic variables
# constants
m, g, d, delta, r, R, L_0, L_1, alpha, c, k, b, phi = sym.symbols('m, g, d, delta, r, R, L_0, L_1, alpha, c, k, b, phi')

# system variables
x1_eq, x2_eq, x3_eq, V_e = sym.symbols('x1_eq, x2_eq, x3_eq, V_e')

# Declare determined values of a-d
A_1 = (5 / (7 * m)) * ((2 * c * (x3_eq ** 2) / (delta - x1_eq) ** 3) - k)
A_2 = -(5 * b) / (7 * m)
A_3 = (10 * c * x3_eq) / (7 * m * ((delta - x1_eq) ** 2))
B_1 = 1 / (L_0 + (L_1 * np.exp(-alpha * (delta - x1_eq))))
B_2 = -(L_1 * alpha * ((-R * x3_eq) + V_e) * np.exp(-alpha * (delta - x1_eq))) / (
            (L_0 + (L_1 * np.exp(-alpha * (delta - x1_eq)))) ** 2)
B_3 = -R / (L_0 + (L_1 * np.exp(-alpha * (delta - x1_eq))))

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
phi_value = 42  # SHOULD THIS BE RADS??
V_e_value = 36.04

x1_eq_value = 0.47861
x2_eq_value = 0
x3_eq_eqn = V_e / R
x3_eq_value = x3_eq_eqn.subs([(V_e, V_e_value), (R, R_value)])

# Substitute values of the constants into the equations for A_1 -> B_3
a_value = a.subs([(M, M_value), (m, m_value)])
b_value = b.subs([(M, M_value), (m, m_value), (g, g_value)])
c_value = c.subs([(M, M_value), (m, m_value), (g, g_value), (ell, ell_value)])
d_value = d.subs([(M, M_value), (m, m_value), (g, g_value), (ell, ell_value)])

# Declare additional symbols from the transfer function
s, t = sym.symbols('s, t')
w = sym.symbols('w', real=True)  # w represents omega
# a-d must be positive for inverse Laplace transform to function correctly
a, b, c, d = sym.symbols('a:d', real=True, positive=True)

# Declare transfer functions G_theta and G_x
# from derivation on additional notes 2
G_theta = -c / (s ** 2 - d)
G_x = ((a * s ** 2) - (a * d) + (b * c)) / (s ** 4 - (s ** 2 * d))

# Perform an impulse (kick), step (push) and frequency (shake) response
# Kick -> Dirac pulse: F_s = 1
# Push -> Step response: F_s = 1/s
# Shake -> Freq. response: F_s = 1/(s**2 + 1)
F_s_kick = 1
F_s_push = 1 / s
F_s_shake = w / (s ** 2 + w ** 2)

# G_theta responses for kick, push and shake
X3_s_kick = G_theta * F_s_kick
x3_t_kick = sym.inverse_laplace_transform(X3_s_kick, s, t)

X3_s_push = G_theta * F_s_push
x3_t_push = sym.inverse_laplace_transform(X3_s_push, s, t)

X3_s_shake = G_theta * F_s_shake
x3_t_shake = sym.inverse_laplace_transform(X3_s_shake, s, t, w)

# G_x responses for kick, push and shake
X1_s_kick = G_x * F_s_kick
x1_t_kick = sym.inverse_laplace_transform(X1_s_kick, s, t)

X1_s_push = G_x * F_s_push
x1_t_push = sym.inverse_laplace_transform(X1_s_push, s, t)

X1_s_shake = G_x * F_s_shake
x1_t_shake = sym.inverse_laplace_transform(X1_s_shake, s, t, w)

# Print equations inc. LaTex code
print('x3 response - kick')
sym.pprint(x3_t_kick.simplify())
print(sym.latex(x3_t_kick.simplify()))
print('x3 response - push')
sym.pprint(x3_t_push.simplify())
print(sym.latex(x3_t_push.simplify()))
print('x3 response - shake')
sym.pprint(x3_t_shake.simplify())
print(sym.latex(x3_t_shake.simplify()))

print('x1 response - kick')
sym.pprint(x1_t_kick.simplify())
print(sym.latex(x1_t_kick.simplify()))
print('x1 response - push')
sym.pprint(x1_t_push.simplify())
print(sym.latex(x1_t_push.simplify()))
print('x1 response - shake')
sym.pprint(x1_t_shake.simplify())
print(sym.latex(x1_t_shake.simplify()))
