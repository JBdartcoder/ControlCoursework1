from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


class Car:

    def __init__(self,
                 radius=0.125,                  # r
                 resistance=53,                 # R
                 nominal_inductance=0.120,      # L_0
                 inductance_constant=0.25,      # L_1
                 inductance_exp_constant=1.2,   # alpha
                 electromagnet_constant=6.815,  # c/1000 so it is in terms of kg
                 spring_stiffness=1880,         # k
                 damper_coeff=10.4,             # b
                 phi=42,                        # phi angle

                 length=2.3,        # set car length = 2.3m as per requirement
                 velocity=5,        # set velocity to 5ms^-1
                 x_position=0.,
                 y_position=0.3,
                 pose=np.deg2rad(5)):          # pose of 5 deg converted to rads
        self.__length = length
        self. __velocity = velocity
        self.__x_position = x_position
        self.__y_position = y_position
        self.__pose_initial = pose

        self.__radius = radius
        self.__resistance = resistance
        self.__nominal_inductance = nominal_inductance
        self.__inductance_constant = inductance_constant
        self.__inductance_exp_constant = inductance_exp_constant
        self.__electromagnet_constant = electromagnet_constant
        self.__spring_stiffness = spring_stiffness
        self.__damper_coeff = damper_coeff
        self.__phi = phi


    def move(self, steering_angle, dt):         # where dt denotes time interval

        # STEP 1: Define the system dynamics
        def system_dynamics(_t,z):
            # Define the system dynamics given in eqns. 2.1a-c
            theta = z[2]    # theta is third part of the function (x,y,theta)
            return [self.__velocity * np.cos(theta),
                    self.__velocity * np.sin(theta),
                    self.__velocity * np.tan(steering_angle) / self.__length]

        # STEP 2: Define the initial conditions, z(0)
        z_initial = [self.__x_position,
                     self.__y_position,
                     self.__pose_initial]

        # STEP 3: Call solve_ivp ~ solves equations
        __num_points = 1000     # Resolution of the solution
        solution = solve_ivp(system_dynamics,
                             [0, dt],
                             z_initial,
                             t_eval=np.linspace(0, dt, __num_points))

        self.__x_position = solution.y[0][-1]   # [-1] represents last element in the array
        self.__y_position = solution.y[1][-1]
        self.__pose_initial = solution.y[2][-1]

        # STEP 4: Plot solutions - plots open one after the other

        # plot x vs time
        plt.plot(solution.t, solution.y[0].T)
        plt.grid()
        plt.xlabel('Time (s)')
        plt.ylabel('x (m)')
        plt.show()

        # plot y vs time
        plt.plot(solution.t, solution.y[1].T)
        plt.grid()
        plt.xlabel('Time (s)')
        plt.ylabel('y (m)')
        plt.show()

        # plot theta vs time
        plt.plot(solution.t, solution.y[2].T)
        plt.grid()
        plt.xlabel('Time (s)')
        plt.ylabel('\u03F4 (radians)')
        plt.show()


# Declare car with steering angle -2 deg for time 2s
# Negative angle as it is clockwise
# Frame of reference determines anticlockwise as +ve
ronald = Car()
ronald.move(np.deg2rad(-2), 2)   # 2 deg converted to rads