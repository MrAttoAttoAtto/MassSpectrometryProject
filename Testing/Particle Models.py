from matplotlib import pyplot as plt
import numpy as np
import math
import time


def show_path(B: float, r: float, V: float, q: float, m: float, timestep: float = 10 ** -9) -> None:
    """
    Shows a matplotlib plot of the path of a particle (mass m, charge q) through a magnetic field of strength B,
    radius (of action) r which is perpendicular to the plane of the page, having been accelerated by a voltage V.

    Why not use vectors/numpy arrays? Because that ends up being 21x slower...

    :param B: The strength of the magnetic field in tesla
    :param r: The radius of the effective area of the magnetic field in metres
    :param V: The accelerating voltage on the particle (before it enters the field) in Volts (assuming start at rest)
    :param q: The charge on the particle in electron charges
    :param m: The mass of the particle in AMUs
    :param timestep: The length of the timesteps that are taken to approximate a constantly changing force
    :return: None
    """

    r_squared = r ** 2

    mass = m * 1.66054e-27
    charge = q * 1.60217662e-19

    # Note that this is constant, as the magnetic field can do no work on the particle as it always contributes
    # a force perpendicular to the direction of motion
    energy = V * charge

    # energy = 1/2*mass*v^2 ==> v = root(2energy/mass)
    speed = math.sqrt(2 * energy / mass)

    # [x, y, z]
    vel_x = 0
    vel_y = speed

    pos_x = 0
    pos_y = -r * 1.5

    points = []
    while pos_x <= 2 * r and pos_y <= 2 * r:
        points.append([pos_x, pos_y])

        # If it is inside the field...
        if pos_x ** 2 + pos_y ** 2 <= r_squared:
            # Lorentz force w/o electric field: cross product
            force_x = charge * vel_y * B
            force_y = charge * -vel_x * B
        else:
            force_x = 0
            force_y = 0

        # F = ma ==> a = F/m
        accel_x = force_x / mass
        accel_y = force_y / mass

        # s = ut + 1/2 * at^2
        disp_x = vel_x * timestep + 1 / 2 * accel_x * timestep ** 2
        pos_x += disp_x
        disp_y = vel_y * timestep + 1 / 2 * accel_y * timestep ** 2
        pos_y += disp_y

        # a = (v - u)/t ==> v = at + u
        vel_x += accel_x * timestep
        vel_y += accel_y * timestep

    fig, ax = plt.subplots()

    ax.plot(*zip(*points))
    circle = plt.Circle((0, 0), r, color="black", fill=False)
    ax.add_artist(circle)
    ax.set(xlim=(-2.2 * r, 2.2 * r), ylim=(-2.2 * r, 2.2 * r))
    fig.show()

    print(points[-1])


def show_path_diff(B: float, r: float, V: float, q: float, m: float, timestep: float = 10 ** -9) -> None:

    mass = m * 1.66054e-27
    charge = q * 1.60217662e-19

    # Note that this is constant, as the magnetic field can do no work on the particle as it always contributes
    # a force perpendicular to the direction of motion
    energy = V * charge

    # energy = 1/2*mass*v^2 ==> v = root(2energy/mass)
    speed = math.sqrt(2 * energy / mass)

    k = charge * B / mass
    t = math.pi/(2*k)

    values = list(np.arange(0, t, timestep))
    points = [[speed*(1-math.cos(k*i))/k for i in values], [speed*math.sin(k*i)/k-r for i in values]]

    fig, ax = plt.subplots()

    ax.plot(points[0], points[1])
    circle = plt.Circle((0, 0), r, color="black", fill=False)
    ax.add_artist(circle)
    ax.set(xlim=(-2.2 * r, 2.2 * r), ylim=(-2.2 * r, 2.2 * r))
    fig.show()


def is_detected(B: float, r: float, V: float, q: float, m: float, back_length: float, slit_width: float) -> bool:

    mass = m * 1.66054e-27
    charge = q * 1.60217662e-19

    # Note that this is constant, as the magnetic field can do no work on the particle as it always contributes
    # a force perpendicular to the direction of motion
    energy = V * charge

    # energy = 1/2*mass*v^2 ==> v = root(2energy/mass)
    speed = math.sqrt(2 * energy / mass)

    # TODO change
    v_x = 0
    v_y = speed

    entry_x = 0
    entry_y = -r

    tangent_gradient = np.inf if v_x == 0 else v_y/v_x

    radius_length = mass * speed / (charge * B)
    radius_gradient = -1/tangent_gradient

    centre_x = radius_length*math.sqrt(1/(1+radius_gradient**2)) + entry_x
    centre_y = radius_gradient * centre_x + entry_y

    a = centre_x
    b = centre_y
    p = radius_length
    k = r ** 2 + a ** 2 + b ** 2 - p ** 2
    exit_y = (b*k + math.sqrt(4 * a**4 * r**2 + 4 * a**2 * b**2 * r**2 - a**2 * k**2)) / (2*(a**2 + b**2))
    exit_x = math.sqrt(r**2 - exit_y**2)

    radius_gradient_2 = np.inf if exit_x == centre_x else (exit_y-centre_y) / (exit_x-centre_x)
    tangent_gradient_2 = -1/radius_gradient_2

    #exit_v_x = speed * math.sqrt(1/(1+tangent_gradient_2**2))
    #exit_v_y = tangent_gradient_2 * exit_v_x

    c = exit_y - tangent_gradient_2 * exit_x
    y = tangent_gradient_2 * (back_length+r) + c

    print(y)

    return abs(y) <= slit_width/2

'''
def show_path_np(B: float, r: float, V: float, e: float, m: float, timestep: float = 10 ** -9) -> None:
    """
    Shows a matplotlib plot of the path of a particle (mass m, charge e) through a magnetic field of strength B,
    radius (of action) r which is perpendicular to the plane of the page, having been accelerated by a voltage V.

    :param B: The strength of the magnetic field in tesla
    :param r: The radius of the effective area of the magnetic field in metres
    :param V: The accelerating voltage on the particle (before it enters the field) in Volts (assuming start at rest)
    :param e: The charge on the particle in electron charges
    :param m: The mass of the particle in AMUs
    :return: None
    """

    r_squared = r ** 2

    mass = m * 1.66054e-27
    charge = e * 1.60217662e-19

    field = np.array([0, 0, B])

    # Note that this is constant, as the magnetic field can do no work on the particle as it always contributes
    # a force perpendicular to the direction of motion
    energy = V * charge

    # energy = 1/2*mass*v^2 ==> v = root(2energy/mass)
    speed = math.sqrt(2 * energy / mass)

    # [x, y, z]
    velocity = np.array([0, speed, 0])
    position = np.array([0, -r, 0])

    points = []

    while max(abs(position)) <= 2 * r:
        points.append([position[0], position[1]])

        # If it is inside the field...
        if position[0] ** 2 + position[1] ** 2 <= r_squared:
            # Lorentz force w/o electric field
            force = charge * np.cross(velocity, field)
        else:
            force = np.array([0, 0, 0])

        # F = ma ==> a = F/m
        acceleration = force / mass

        # s = ut + 1/2 * at^2
        displacement = velocity * timestep + 1 / 2 * acceleration * timestep ** 2
        position += displacement

        # a = (v - u)/t ==> v = at + u
        velocity += acceleration * timestep

    #fig, ax = plt.subplots()

    #ax.plot(*zip(*points))
    #circle = plt.Circle((0, 0), r, color="black", fill=False)
    #ax.add_artist(circle)
    #fig.show()
'''


show_path_diff(0.26, 0.05, 550, 1, 15, 1e-10)
print(is_detected(0.26, 0.05, 550, 1, 15, 0.04, 6.35e-5))
