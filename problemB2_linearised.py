from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym

# Define all necessary symbols to calculate derivatives
x1, x2, x3, x1e, x2e, x3e, c, m, g, k, d, b, V, Ve, R, L0, L1, delta, alpha = sym.symbols(
    'x1, x2, x3, x1e, x2e, x3e, c, m, g, k, d, b, V, Ve, R, L0, L1, delta, alpha')

# Define phi
phi = (c * x3 ** 2 / (delta - x1) ** 2) - (c * x3e ** 2 / (delta - x1e) ** 2) - k * (x1 - x1e) - b * (x2 - x2e)
phi *= 5 / (7 * m)

# Define psi
psi = (V - x3 * R) / (L0 + L1 * sym.exp(-alpha * (delta - x1)))
psi -= (Ve - x3e * R) / (L0 + L1 * sym.exp(-alpha * (delta - x1e)))


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
    This function evaluates the equation
    at the x2=0
    """
    return f.subs([(x2, 0)])


# Evaluate phi at equilibrium point
d_phi_x1_eq = evaluate_at_equilibrium_phi(d_phi_x1)
d_phi_x2_eq = evaluate_at_equilibrium_phi(d_phi_x2)
d_phi_x3_eq = evaluate_at_equilibrium_phi(d_phi_x3)

# End of derivatives calculation


# Constants defined for linearised system
A_1 = d_phi_x1_eq
A_2 = d_phi_x2_eq
A_3 = d_phi_x3_eq
B_1 = d_psi_V
B_2 = d_psi_x1
B_3 = d_psi_x3


class WoodenBall:

    # Declared the values of constants
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
                 x_1_value=0.46,
                 x_2_value=0,
                 x_3_value=0,
                 x_1_e_value=0.47861,  # x_1
                 x_2_e_value=0,  # x_2
                 x_3_e_value=0.68,  # x_3
                 V_e_value=36.04
                 ):
        self.__mass = m_value
        self.__gravity = g_value
        self.__natural_length = d_value
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
            x_1a = z[0]
            x_2a = z[1]
            x_3a = z[2]
            V_value = voltage

            # Subbing the values into the constants: A_1, A_2 ... etc
            def evaluate_constants(f):
                """
                This function subs in
                the constant values into
                A1, A2, A3, B1, B2 and B3
                """
                return f.subs([(m, self.__mass), (R, self.__resistance), (L0, self.__normal_inductance),
                               (L1, self.__inductance_constant), (alpha, self.__inductance_exp_constant),
                               (c, self.__electromagnet_constant),
                               (delta, self.__magnet_position), (k, self.__spring_stiffness),
                               (b, self.__damper_coeff), (x3, x_3a), (x1, x_1a), (V, V_value)])

            A_1_sub = evaluate_constants(A_1)
            A_2_sub = evaluate_constants(A_2)
            A_3_sub = evaluate_constants(A_3)
            B_1_sub = evaluate_constants(B_1)
            B_2_sub = evaluate_constants(B_2)
            B_3_sub = evaluate_constants(B_3)

            # Declaring the bar variables
            x_1_bar = x_1a - self.x_1_e
            x_2_bar = x_2a - self.x_2_e
            x_3_bar = x_3a - self.x_3_e
            V_bar = voltage - self.V_e

            # The equations are calculated using the variables and returned
            return [x_2_bar, (A_1_sub * x_1_bar) + (A_2_sub * x_2_bar) + (A_3_sub * x_3_bar),
                    (B_1_sub * V_bar) + (B_2_sub * x_1_bar) + (B_3_sub * x_3_bar)]

        # STEP 2: Define the initial conditions, z(0)
        z_initial = [self.x_1,
                     self.x_2,
                     self.x_3]

        # STEP 3: Call solve_ivp ~ solves equations
        __num_points = 1000  # Resolution of the solution
        solution = solve_ivp(system_dynamics,
                             [0, dt],
                             z_initial,
                             t_eval=np.linspace(0, dt, __num_points))

        self.x_1 = solution.y[0][-1]
        self.x_2 = solution.y[1][-1]
        self.x_3 = solution.y[2][-1]

        # STEP 4: Plot solutions

        # Subplot is used to plot the 3 graphs in one figure
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

        ax1.plot(solution.t, solution.y[0])
        ax1.grid()
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('$x_1$ (m)')

        ax2.plot(solution.t, solution.y[1])
        ax2.grid()
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('$x_2$ (m/s)')

        ax3.plot(solution.t, solution.y[2])
        ax3.grid()
        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('$x_3$ (A)')

        plt.show()


woody = WoodenBall()
woody.move(36.04, 2)
