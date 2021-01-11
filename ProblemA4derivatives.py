import sympy as sym

print('This code is used to produces the derivatives for phi and psi which is used to answer Problem A4.')
print()

# Define all necessary symbols
x1, x2, x3, x1e, x2e, x3e, c, m, g, k, d, b, theta, V, Ve, R, L0, L1, delta, alpha = sym.symbols(
    'x1, x2, x3, x1e, x2e, x3e, c, m, g, k, d, b, theta, V, Ve, R, L0, L1, delta, alpha')

# Define phi
phi = (c * x3 ** 2 / (delta - x1) ** 2) - (c*x3e**2 / (delta-x1e)**2) - k * (x1-x1e) - b * (x2-x2e)
phi *= 5 / (7*m)

# Printing out phi
print('Phi = ')
sym.pprint(phi)
print()

# Define psi
psi = (V - x3 * R) / (L0 + L1 * sym.exp(-alpha * (delta - x1)))
psi -= (Ve - x3e * R) / (L0 + L1 * sym.exp(-alpha * (delta - x1e)))

# Printing put psi
print('Psi = ')
sym.pprint(psi)
print()

# Differentiate phi wrt x1, x2, x3
d_phi_x1 = phi.diff(x1)
d_phi_x2 = phi.diff(x2)
d_phi_x3 = phi.diff(x3)


# Differentiate psi wrt V, x1, x3
d_psi_V = psi.diff(V)
d_psi_x1 = psi.diff(x1)
d_psi_x3 = psi.diff(x3)


def evaluate_at_equilibrium_phi(f):
    """
    This function evaluates the phi
    at x2=0
    """
    return f.subs([(x2, 0)])


# Evaluate phi at equilibrium point
d_phi_x1_eq = evaluate_at_equilibrium_phi(d_phi_x1)
d_phi_x2_eq = evaluate_at_equilibrium_phi(d_phi_x2)
d_phi_x3_eq = evaluate_at_equilibrium_phi(d_phi_x3)

# Psi does not need to be evaluated at equilibrium point

# Print partial derivatives of phi
print('d_phi_x1=')
sym.pprint(d_phi_x1_eq)
print()
print('d_phi_x2=')
sym.pprint(d_phi_x2_eq)
print()
print('d_phi_x3=')
sym.pprint(d_phi_x3_eq)
print()

# Print partial derivatives of psi
print('d_psi_V=')
sym.pprint(d_psi_V)
print()
print('d_psi_x1=')
sym.pprint(d_psi_x1)
print()
print('d_psi_x3=')
sym.pprint(d_psi_x3)
