import sympy as sym
import numpy as np

# Define all necessary symbols
x1, x2, x3, x4, x1e, x2e, x3e, x4e, c, m, g, k, d, b, theta, V, Ve, R, L0, L1, sig, alp = sym.symbols(
    'x1, x2, x3, x4, x1e, x2e, x3e, x4e, c, m, g, k, d, b, theta, V, Ve R, L0, L1, sig, alp')

X1, V, p, u, s, t, sig1, d, q = sym.symbols('X1, V, p, u, s, t, sig1, d, q')
G = p*u
G /= ((s**2 - (t*s) - sig1)*(s**2-d)) - p*q
sym.pprint(G)
print()

print(sym.latex(G.factor()))
print()

"""
x1. = x2
x2. = phi
x3. = x4
x4. = psi
"""
# Define phi
phi = (c * x3 ** 2 / (sig - x1) ** 2) - (c*x3e**2 / (sig-x1e)**2) - k * (x1-x1e) - b * (x2-x2e)
phi *= 5 / (7*m)
# sym.pprint(phi)

# Define psi
psi = (V - x3 * R) / (L0 + L1 * sym.exp(-alp * (sig - x1)))
psi -= (Ve - x3e * R) / (L0 + L1 * sym.exp(-alp * (sig - x1e)))
sym.pprint(psi)

# Differentiate phi wrt x1, x2, x3
d_phi_x1 = phi.diff(x1)
d_phi_x2 = phi.diff(x2)
d_phi_x3 = phi.diff(x3)

print()

# Differentiate psi wrt V, x1, x3
d_psi_V = psi.diff(V)
d_psi_x1 = psi.diff(x1)
d_psi_x3 = psi.diff(x3)


def evaluate_at_equilibrium_phi(f):
    """
    This function evaluates the equation
    at the x1=0, x2=0 and x3=0
    """
    return f.subs([(x2, 0)])


def evaluate_at_equilibrium_psi(f):
    """
    This function evaluates the equation
    at the V=0, x1=0 and x3=0
    """
    return f.subs()


# Evaluate phi at equilibrium point
d_phi_x1_eq = evaluate_at_equilibrium_phi(d_phi_x1)
d_phi_x2_eq = evaluate_at_equilibrium_phi(d_phi_x2)
d_phi_x3_eq = evaluate_at_equilibrium_phi(d_phi_x3)

# Evaluate psi at equilibrium point
# d_psi_V_eq = evaluate_at_equilibrium_psi(d_psi_V)
# d_psi_x1_eq = evaluate_at_equilibrium_psi(d_psi_x1)
# d_psi_x3_eq = evaluate_at_equilibrium_psi(d_psi_x3)


# Print partial derivatives of phi
print('d_phi_x1:')
sym.pprint(d_phi_x1_eq)
print('d_phi_x2:')
sym.pprint(d_phi_x2_eq)
print('d_phi_x3:')
sym.pprint(d_phi_x3_eq)


# Print partial derivatives of psi
print('d_psi_V:')
sym.pprint(d_psi_V)
print('d_psi_x1:')
sym.pprint(d_psi_x1)
print('d_psi_x2:')
sym.pprint(d_psi_x3)

# End
# Simon added this code
# John is commenting to test Git
