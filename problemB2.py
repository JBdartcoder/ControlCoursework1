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
                 inductance_constant=0.25,      # L_1
                 alpha=1.2,   # alpha
                 c=6.815,  # c/1000 so it is in terms of kg
                 k=1880,         # k
                 b=10.4,             # b
                 phi=42,                        # phi angle
                 # v=36.04,                 # V
                 x_1_position=0.45,        # x_1
                 x_2=0,                 # x_2
                 x_3=0                  # x_3
                 ):

        self.__mass = m
        self.__gravity = g
        self.__natural_lenth = d
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
        # self.__voltage = v
        self.x_1 = x_1_position
        self.x_2 = x_2
        self.x_3 = x_3

    def move(self, voltage, dt):         # where dt denotes time interval

        # STEP 1: Define the system dynamics
        def system_dynamics(_t, z):
            x_2 = z[1]
            x_3 = z[2]
            return [x_2,
                    (5 / 7 * self.__mass) * ((self.__electromagnet_constant * (x_3 ** 2) / (self.__magnet_position - self.x_1) ** 2)
                                             + (self.__mass * self.__gravity * np.sin(self.__phi)) - (self.__damper_coeff * x_2) -
                                             (self.__spring_stiffness*(self.x_1 - self.__natural_lenth))),
                    (voltage - (x_3 * self.__resistance)) /
                    (self.__nominal_inductance + (self.__inductance_constant * np.exp(-self.__inductance_exp_constant*(self.__magnet_position - self.x_1))))
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
ronald.move(36.04, 2)   # 2 deg converted to rads
