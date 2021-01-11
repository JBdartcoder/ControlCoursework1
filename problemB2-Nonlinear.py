from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


class WoodenBall:

    def __init__(self,
                 m=0.425,                    # m
                 g=9.81,                  # g
                 d=0.42,           # d
                 delta=0.65,          # delta
                 r=0.125,                  # r
                 R=53,                 # R
                 nominal_inductance=0.120,      # L_0
                 inductance_constant=0.025,      # L_1
                 alpha=1.2,   # alpha
                 c=6.815,  # c/1000 so it is in terms of kg
                 k=1880,         # k
                 b=10.4,             # b
                 phi=0.73303,                        # phi angle
                 x_1=0.46,        # x_1
                 x_2=0,                 # x_2
                 x_3=0                  # x_3
                 ):

        self.__mass = m
        self.__gravity = g
        self.__natural_length = d
        self.__magnet_position = delta
        self.__radius = r
        self.__resistance = R
        self.__nominal_inductance = nominal_inductance
        self.__inductance_constant = inductance_constant
        self.__inductance_exp_constant = alpha
        self.__electromagnet_constant = c
        self.__spring_stiffness = k
        self.__damper_coeff = b
        self.__phi = phi
        self.x_1 = x_1
        self.x_2 = x_2
        self.x_3 = x_3

    def move(self, voltage, dt):         # where dt denotes time interval

        # STEP 1: Define the system dynamics
        def system_dynamics(_t, z):
            x_1a = z[0]
            x_2a = z[1]
            x_3a = z[2]
            return [x_2a,
                    (5 / 7 * self.__mass) * ((self.__electromagnet_constant * (x_3a ** 2) / (self.__magnet_position - x_1a) ** 2)
                                             + (self.__mass * self.__gravity * np.sin(self.__phi)) - (self.__damper_coeff * x_2a) -
                                             (self.__spring_stiffness*(x_1a - self.__natural_length))),
                    (voltage - (x_3a * self.__resistance)) /
                    (self.__nominal_inductance + (self.__inductance_constant * np.exp(-self.__inductance_exp_constant*(self.__magnet_position - x_1a))))
                    ]

        # STEP 2: Define the initial conditions, z(0)
        z_initial = [self.x_1,
                     self.x_2,
                     self.x_3]

        # STEP 3: Call solve_ivp ~ solves equations
        __num_points = 1000     # Resolution of the solution
        solution = solve_ivp(system_dynamics,
                             [0, dt],
                             z_initial,
                             t_eval=np.linspace(0, dt, __num_points))

        self.x_1 = solution.y[0][-1]   # [-1] represents last element in the array
        self.x_2 = solution.y[1][-1]
        self.x_3 = solution.y[2][-1]

        # STEP 4: Plot solutions - plots open one after the other

        """
        # plot x vs time
        plt.plot(solution.t, solution.y[0].T)
        plt.grid()
        plt.xlabel('Time (s)')
        plt.ylabel('$x_1$ (m)')
        plt.show()

        # plot y vs time
        plt.plot(solution.t, solution.y[1].T)
        plt.grid()
        plt.xlabel('Time (s)')
        plt.ylabel('$x_2$ (m/s)')
        plt.show()

        # plot theta vs time
        plt.plot(solution.t, solution.y[2].T)
        plt.grid()
        plt.xlabel('Time (s)')
        plt.ylabel('$x_3$ (A)')
        plt.show()
        """

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


ronald = WoodenBall()
ronald.move(36.04, 2)
