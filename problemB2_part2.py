from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym

# Define all necessary symbols
x1, x2, x3, x1e, x2e, x3e, c, m, g, k, d, b, V, Ve, R, L0, L1, delta, alpha = sym.symbols(
    'x1, x2, x3, x1e, x2e, x3e, c, m, g, k, d, b, V, Ve, R, L0, L1, delta, alpha')

X1, V, p, u, s, t, sig1, d, q = sym.symbols('X1, V, p, u, s, t, sig1, d, q')
G = p * u
G /= ((s ** 2 - (t * s) - sig1) * (s ** 2 - d)) - p * q
# sym.pprint(G)
print()

# print(sym.latex(G.factor()))
print()

"""
x1. = x2
x2. = phi
x3. = psi
"""
# Define phi
phi = (c * x3 ** 2 / (delta - x1) ** 2) - (c * x3e ** 2 / (delta - x1e) ** 2) - k * (x1 - x1e) - b * (x2 - x2e)
phi *= 5 / (7 * m)
# sym.pprint(phi)

# Define psi
psi = (V - x3 * R) / (L0 + L1 * sym.exp(-alpha * (delta - x1)))
psi -= (Ve - x3e * R) / (L0 + L1 * sym.exp(-alpha * (delta - x1e)))
# sym.pprint(psi)

# Differentiate phi wrt x1, x2, x3
d_phi_x1 = phi.diff(x1)
d_phi_x2 = phi.diff(x2)
d_phi_x3 = phi.diff(x3)

# print()

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
    return f.subs([])


# Evaluate phi at equilibrium point
d_phi_x1_eq = evaluate_at_equilibrium_phi(d_phi_x1)
d_phi_x2_eq = evaluate_at_equilibrium_phi(d_phi_x2)
d_phi_x3_eq = evaluate_at_equilibrium_phi(d_phi_x3)

# Evaluate psi at equilibrium point
d_psi_V_eq = evaluate_at_equilibrium_psi(d_psi_V)
d_psi_x1_eq = evaluate_at_equilibrium_psi(d_psi_x1)
d_psi_x3_eq = evaluate_at_equilibrium_psi(d_psi_x3)

"""
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
print('d_psi_x3:')
sym.pprint(d_psi_x3)
"""

# End of derivatives calculation


# Constants defined for linearised system

A_1 = d_phi_x1_eq
A_2 = d_phi_x2_eq
A_3 = d_phi_x3_eq
B_1 = d_psi_V
B_2 = d_psi_x1
B_3 = d_psi_x3


# Define linearised system variables
# x_1_bar, x_2_bar, x_3_bar, V_bar = sym.symbols('x_1_bar, x_2_bar, x_3_bar, V_bar')


class WoodenBall:

    def __init__(self,
                 m_value=0.425,  # m
                 g_value=9.81,  # g
                 d_value=0.42,  # d
                 delta_value=0.65,  # delta
                 r_value=0.125,  # r
                 R_value=53,  # R
                 normal_inductance_value=0.120,  # L_0
                 inductance_constant_value=0.025,  # L_1
                 alpha_value=1.2,  # alpha
                 c_value=6.815,  # c/1000 so it is in terms of kg
                 k_value=1880,  # k
                 b_value=10.4,  # b
                 phi_value=0.73303,  # phi angle
                 # v=36.04,                 # V
                 x_1_value=0,
                 x_2_value=0,
                 x_3_value=0,
                 x_1_e_value=0.47861,  # x_1
                 x_2_e_value=0,  # x_2
                 x_3_e_value=0,  # x_3
                 V_e_value=36.04
                 ):
        self.__mass = m_value
        self.__gravity = g_value
        self.__natural_lenth = d_value
        self.__magnet_position = delta_value
        self.__radius = r_value
        self.__resistance = R_value
        self.__normal_inductance = normal_inductance_value
        self.__inductance_constant = inductance_constant_value
        self.__inductance_exp_constant = alpha_value
        self.__electromagnet_constant = c_value
        self.__spring_stiffness = k_value
        self.__damper_coeff = b_value
        self.__phi = phi_value
        self.x_1 = x_1_value
        self.x_2 = x_2_value
        self.x_3 = x_3_value
        self.x_1_e = x_1_e_value
        self.x_2_e = x_2_e_value
        self.x_3_e = x_3_e_value
        self.V_e = V_e_value

    def move(self, voltage, dt):  # where dt denotes time interval

        # STEP 1: Define the system dynamics
        def system_dynamics(_t, z):
            x_1 = z[0]
            x_2 = z[1]
            x_3 = z[2]
            V_value = voltage

            # print(x_1)
            # print(x_3)

            # Subbing the values into the constants: A_1, A_2 ... etc
            A_1_sub = A_1.subs([(m, self.__mass), (c, self.__electromagnet_constant), (delta, self.__magnet_position),
                                (k, self.__spring_stiffness), (x3, x_3), (x1, x_1)])
            A_2_sub = A_2.subs([(b, self.__damper_coeff), (m, self.__mass)])
            A_3_sub = A_3.subs([(c, self.__electromagnet_constant), (m, self.__mass), (delta, self.__magnet_position),
                                (x3, x_3), (x1, x_1)])
            B_1_sub = B_1.subs([(L0, self.__normal_inductance), (L1, self.__inductance_constant),
                                (alpha, self.__inductance_exp_constant), (delta, self.__magnet_position), (x1, x_1)])
            B_2_sub = B_2.subs([(L1, self.__inductance_constant), (alpha, self.__inductance_exp_constant),
                                (R, self.__resistance), (delta, self.__magnet_position),
                                (L0, self.__normal_inductance), (x3, x_3), (x1, x_1), (V, V_value)])
            B_3_sub = B_3.subs([(R, self.__resistance), (L0, self.__normal_inductance),
                                (L1, self.__inductance_constant), (alpha, self.__inductance_exp_constant),
                                (delta, self.__magnet_position), (x1, x_1)])

            """
            print(A_1_sub)
            print(A_2_sub)
            print(A_3_sub)
            print(B_1_sub)
            print(B_2_sub)
            print(B_3_sub)


            x_1_bar_dot = x_2_bar
            x_2_bar_dot = (A_1 * x_1_bar) - (A_2 * x_2_bar) + (A_3 * x_3_bar)
            x_3_bar_dot = (B_1 * V_bar) - (B_2 * x_1_bar) - (B_3 * x_3_bar)
            """

            x_1_bar = self.x_1 - self.x_1_e
            x_2_bar = self.x_2 - self.x_2_e
            x_3_bar = self.x_3 - self.x_3_e
            V_bar = voltage - self.V_e
            return [x_2_bar, (A_1_sub * x_1_bar) - (A_2_sub * x_2_bar) + (A_3_sub * x_3_bar),
                    (B_1_sub * V_bar) - (B_2_sub * x_1_bar) - (B_3_sub * x_3_bar)]

        # STEP 2: Define the initial conditions, z(0)
        z_initial = [0.47861,
                     self.x_2,
                     self.x_3]

        # STEP 3: Call solve_ivp ~ solves equations
        __num_points = 1000  # Resolution of the solution
        solution = solve_ivp(system_dynamics,
                             [0, dt],
                             z_initial,
                             t_eval=np.linspace(0, dt, __num_points))

        self.x_1 = solution.y[0][-1]  # [-1] represents last element in the array
        self.x_2 = solution.y[1][-1]
        self.x_3 = solution.y[2][-1]

        # STEP 4: Plot solutions - plots open one after the other

        # plot x vs time
        plt.plot(solution.t, solution.y[0].T)
        plt.grid()
        plt.xlabel('Time (s)')
        plt.ylabel('x_1 (m)')
        plt.show()

        # plot y vs time
        plt.plot(solution.t, solution.y[1].T)
        plt.grid()
        plt.xlabel('Time (s)')
        plt.ylabel('x_2')
        plt.show()

        # plot theta vs time
        plt.plot(solution.t, solution.y[2].T)
        plt.grid()
        plt.xlabel('Time (s)')
        plt.ylabel('x_3')
        plt.show()


# Declare car with steering angle -2 deg for time 2s
# Negative angle as it is clockwise
# Frame of reference determines anticlockwise as +ve
ronald = WoodenBall()
ronald.move(36.04, 2)  # 2 deg converted to rads
