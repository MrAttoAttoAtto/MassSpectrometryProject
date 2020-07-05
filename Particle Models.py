import math

import matplotlib.lines
import matplotlib.patches
import numpy as np
from matplotlib import pyplot as plt

import random


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
    t = math.pi / (2 * k)

    values = list(np.arange(0, t, timestep))
    points = [[speed * (1 - math.cos(k * i)) / k for i in values], [speed * math.sin(k * i) / k - r for i in values]]

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

    tangent_gradient = np.inf if v_x == 0 else v_y / v_x

    radius_length = mass * speed / (charge * B)
    radius_gradient = -1 / tangent_gradient

    centre_x = radius_length * math.sqrt(1 / (1 + radius_gradient ** 2)) + entry_x
    centre_y = radius_gradient * centre_x + entry_y

    a = centre_x
    b = centre_y
    p = radius_length
    k = r ** 2 + a ** 2 + b ** 2 - p ** 2
    exit_y = (b * k + math.sqrt(4 * a ** 4 * r ** 2 + 4 * a ** 2 * b ** 2 * r ** 2 - a ** 2 * k ** 2)) / (
                2 * (a ** 2 + b ** 2))
    exit_x = math.sqrt(r ** 2 - exit_y ** 2)

    radius_gradient_2 = np.inf if exit_x == centre_x else (exit_y - centre_y) / (exit_x - centre_x)
    tangent_gradient_2 = -1 / radius_gradient_2

    # exit_v_x = speed * math.sqrt(1/(1+tangent_gradient_2**2))
    # exit_v_y = tangent_gradient_2 * exit_v_x

    c = exit_y - tangent_gradient_2 * exit_x
    y = tangent_gradient_2 * (back_length + r) + c

    print(y)

    return abs(y) <= slit_width / 2


def is_detected_angle(B: float, r: float, V: float, q: float, m: float, entry_length: float, exit_length: float,
                      angle: float, slit_width: float):
    """
    Checks is a given fired particle would be detected by the machine

    :param B:
    :param r:
    :param V:
    :param q:
    :param m:
    :param front_length:
    :param angle:
    :param back_length:
    :param slit_width:
    :return:
    """

    mass = m * 1.66054e-27
    charge = q * 1.60217662e-19

    angle_radians = 2 * math.pi * (90 - angle) / 360

    # Note that this is constant, as the magnetic field can do no work on the particle as it always contributes
    # a force perpendicular to the direction of motion
    energy = V * charge

    # energy = 1/2*mass*v^2 ==> v = root(2energy/mass)
    speed = math.sqrt(2 * energy / mass)

    path_gradient = np.inf if angle == 0 else math.tan(angle_radians)

    m = path_gradient
    c = -(entry_length + r)
    op = (lambda x, y: x - y) if angle >= 0 else (lambda x, y: x + y)
    entry_x = (op(-2 * c * m, math.sqrt(4 * c ** 2 * m ** 2 - 4 * (m ** 2 + 1) * (c ** 2 - r ** 2)))) / (
            2 * (m ** 2 + 1)) if path_gradient != np.inf else 0
    entry_y = -math.sqrt(r ** 2 - entry_x ** 2)

    tangent_gradient = path_gradient

    radius_length = mass * speed / (charge * B)
    radius_gradient = -1 / tangent_gradient

    centre_x = radius_length * math.sqrt(1 / (1 + radius_gradient ** 2)) + entry_x
    centre_y = radius_gradient * centre_x + (entry_y - radius_gradient * entry_x)

    a = centre_x
    b = centre_y
    p = radius_length
    k = r ** 2 + a ** 2 + b ** 2 - p ** 2
    determinant = math.sqrt(4 * a ** 4 * r ** 2 + 4 * a ** 2 * b ** 2 * r ** 2 - a ** 2 * k ** 2)

    exit_y = (b * k + determinant) / (2 * (a ** 2 + b ** 2))
    if abs(exit_y - entry_y) < 1e-12:
        # It is the other solution!
        exit_y = (b * k - determinant) / (2 * (a ** 2 + b ** 2))

    exit_x = math.sqrt(r ** 2 - exit_y ** 2)
    if abs(((exit_x - centre_x) ** 2 + (exit_y - centre_y) ** 2) - radius_length ** 2) >= 1e-12 or \
            abs(exit_x - entry_x) < 1e-12:
        # It is the other x solution!
        exit_x = -exit_x

    radius_gradient_2 = np.inf if exit_x == centre_x else (exit_y - centre_y) / (exit_x - centre_x)
    tangent_gradient_2 = -1 / radius_gradient_2

    # exit_v_x = speed * math.sqrt(1/(1+tangent_gradient_2**2))
    # exit_v_y = tangent_gradient_2 * exit_v_x

    # x^2 + y^2 = distance^2 (y = mx + c simplified)
    distance_right = (exit_x + 1) ** 2 + (tangent_gradient_2 + exit_y) ** 2
    distance_left = (exit_x - 1) ** 2 + (exit_y - tangent_gradient_2) ** 2

    going_right = distance_right >= distance_left

    if not going_right:
        return False

    c = exit_y - tangent_gradient_2 * exit_x
    final_y = tangent_gradient_2 * (exit_length + r) + c

    return abs(final_y) <= slit_width / 2


def show_path_diff_angle(B: float, r: float, V: float, q: float, m: float, entry_length: float, exit_length: float,
                         angle: float, slit_width: float):
    """

    :param B:
    :param r:
    :param V:
    :param q:
    :param m:
    :param front_length:
    :param back_length:
    :param angle:
    :param slit_width:
    :return:
    """

    mass = m * 1.66054e-27
    charge = q * 1.60217662e-19

    angle_radians = 2 * math.pi * (90 - angle) / 360

    # Note that this is constant, as the magnetic field can do no work on the particle as it always contributes
    # a force perpendicular to the direction of motion
    energy = V * charge

    # energy = 1/2*mass*v^2 ==> v = root(2energy/mass)
    speed = math.sqrt(2 * energy / mass)

    path_gradient = np.inf if angle == 0 else math.tan(angle_radians)

    m = path_gradient
    c = -(entry_length + r)
    op = (lambda x, y: x - y) if angle >= 0 else (lambda x, y: x + y)
    entry_x = (op(-2 * c * m, math.sqrt(4 * c ** 2 * m ** 2 - 4 * (m ** 2 + 1) * (c ** 2 - r ** 2)))) / (
            2 * (m ** 2 + 1)) if path_gradient != np.inf else 0
    entry_y = -math.sqrt(r ** 2 - entry_x ** 2)

    tangent_gradient = path_gradient

    radius_length = mass * speed / (charge * B)
    radius_gradient = -1 / tangent_gradient

    centre_x = radius_length * math.sqrt(1 / (1 + radius_gradient ** 2)) + entry_x
    centre_y = radius_gradient * centre_x + (entry_y - radius_gradient * entry_x)

    a = centre_x
    b = centre_y
    p = radius_length
    k = r ** 2 + a ** 2 + b ** 2 - p ** 2
    determinant = math.sqrt(4 * a ** 4 * r ** 2 + 4 * a ** 2 * b ** 2 * r ** 2 - a ** 2 * k ** 2)

    exit_y = (b * k + determinant) / (2 * (a ** 2 + b ** 2))
    if abs(exit_y - entry_y) < 1e-12:
        # It is the other solution!
        exit_y = (b * k - determinant) / (2 * (a ** 2 + b ** 2))

    exit_x = math.sqrt(r ** 2 - exit_y ** 2)
    if abs(((exit_x - centre_x) ** 2 + (exit_y - centre_y) ** 2) - radius_length ** 2) >= 1e-12 or abs(
            exit_x - entry_x) < 1e-12:
        # It is the other x solution!
        exit_x = -exit_x

    radius_gradient_2 = np.inf if exit_x == centre_x else (exit_y - centre_y) / (exit_x - centre_x)
    tangent_gradient_2 = -1 / radius_gradient_2

    fig, ax = plt.subplots()

    circle = plt.Circle((0, 0), r, color="black", fill=False)
    ax.add_artist(circle)

    exit_above_line = exit_y > (radius_gradient * exit_x + (entry_y - radius_gradient * entry_x))
    op = (lambda x, y: x - y) if exit_above_line else (lambda x, y: x + y)
    end_angle = op((180 - angle), 180 * math.acos(
        1 - ((entry_x - exit_x) ** 2 + (entry_y - exit_y) ** 2) / (2 * radius_length ** 2)) / math.pi)
    arc_path = matplotlib.patches.Arc((centre_x, centre_y), radius_length * 2, radius_length * 2, 0,
                                      end_angle, 180 - angle, color="green")

    if angle != 0:
        x1 = np.linspace(0, entry_x, 1000)
        y1 = tangent_gradient * x1 - r - entry_length
    else:
        x1 = np.zeros(1000)
        y1 = np.linspace(-r - entry_length, -r, 1000)

    # x^2 + y^2 = distance^2 (y = mx + c simplified)
    distance_right = (exit_x + 1) ** 2 + (tangent_gradient_2 + exit_y) ** 2
    distance_left = (exit_x - 1) ** 2 + (exit_y - tangent_gradient_2) ** 2

    going_right = distance_right >= distance_left

    end_x = r + exit_length if going_right else -r - exit_length
    x2 = np.linspace(exit_x, end_x, 1000)
    y2 = tangent_gradient_2 * x2 + (exit_y - tangent_gradient_2 * exit_x)

    ax.plot(x1, y1)
    ax.plot(x2, y2)

    ax.vlines(r + exit_length, -slit_width / 2, -r,  color="purple")
    ax.vlines(r + exit_length, slit_width / 2, r, color="purple")

    ax.add_patch(arc_path)
    lim = max(1.2 * r + entry_length, 1.2 * r + exit_length)
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim))
    fig.show()


show_path_diff_angle(0.3020026904, 0.05, 550, 1, 20, 0.05, 0.05, 0, 0.01)
#cool = [[random.random()+1, random.random()*0.05+0.01] for i in range(10000)]
#for i in range(10000):
#    is_detected_angle(0.6752985451, cool[i][1], 550, 1, 100, 0.05, 0.05, cool[i][0], 6.7e-5)
