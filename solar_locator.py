######################################################################
# This file copyright the Georgia Institute of Technology
#
# Permission is given to students to use or modify this file (only)
# to work on their assignments.
#
# You may NOT publish this file or make it available to others not in
# the course.
#
######################################################################

#These import statements give you access to library functions which you may
# (or may not?) want to use.
from math import *
import random
from body import *
from solar_system import *
import math

def estimate_next_pos(solar_system, gravimeter_measurement, other=None):
    """
    Estimate the next (x,y) position of the satelite.
    This is the function you will have to write for part A.

    :param solar_system: SolarSystem
        A model of the positions, velocities, and masses
        of the planets in the solar system, as well as the sun.
    :param gravimeter_measurement: float
        A floating point number representing
        the measured magnitude of the gravitation pull of all the planets
        felt at the target satellite at that point in time.
    :param other: any
        This is initially None, but if you return an OTHER from
        this function call, it will be passed back to you the next time it is
        called, so that you can use it to keep track of important information
        over time. (We suggest you use a dictionary so that you can store as many
        different named values as you want.)
    :return:
        estimate: Tuple[float, float]. The (x,y) estimate of the target satellite at the next timestep
        other: any. Any additional information you'd like to pass between invocations of this function
        optional_points_to_plot: List[Tuple[float, float, float]].
            A list of tuples like (x,y,h) to plot for the visualization
    """

    # example of how to get the gravity magnitude at a body in the solar system:
    particle = Body(r=[1*AU, 1*AU], v=[0, 0], mass=0, measurement_noise=0)
    particle_gravimeter_measurement = particle.compute_gravity_magnitude(planets=solar_system.planets)

    # You must return a tuple of (x,y) estimate, and OTHER (even if it is NONE)
    # in this order for grading purposes.

    xy_estimate = (0, 0)  # Sample answer, (X,Y) as a tuple.

    # TODO - remove this canned answer which makes this template code
    # pass one test case once you start to write your solution....
    # Step 1 - Initialize - want to generate set of N random particles in space that cover object
    # search space is [+-4, +-4]
    # randomly generate points 20 radii from 0 - 4
    # stagger 20 points on each radii (200 points total)
    samplePoints = []
    timer = other
    if not timer:
        for i in range(20):
            for j in range(20):
                radius = (i+1) * 1/10 # this gives us radial distance from sun
                theta = random.random() * 2 * pi # this gives us angle
                x = radius * cos(theta)
                y = radius * sin(theta)
                thisPoint = Particle((x,y))
                samplePoints.append(thisPoint)
                timer = Timer(0, samplePoints)
    else:
        for i in range(len(timer.particles)):
            samplePoints.append(Particle(timer.particles[i].position))
        # for i in range(10): # adds random component to points
        #     radius = random.random() * 4
        #     theta = random.random() * 2 * pi  # this gives us angle
        #     x = radius * cos(theta)
        #     y = radius * sin(theta)
        #     thisPoint = Particle((x, y))
        #     samplePoints.append(thisPoint)

    # for i in range(len(samplePoints)):
    #     print(samplePoints[i].position)

    # Step 2 - Assign importance weights - how likely particle is to be object - assign based on gravimeter
    # first need gravimeter measure
    objectMeasure = gravimeter_measurement
    #print('objectMeasure:', objectMeasure)
    # need to assign score to each point -- closer to 0, the closer the score
    for point in samplePoints:
        pointBody = Body(r=[AU * point.position[0], AU * point.position[1]], v=[0, 0], mass=0, measurement_noise=0)
        pointMeasure = pointBody.compute_gravity_magnitude(planets=solar_system.planets)
        point.weight = abs(pointMeasure - objectMeasure)
        point.measure = pointMeasure
        #print(pointMeasure)



    # Step 3 - Resample - generate N new points based on importance weights of last set
    # keep the 25 best particles??? -- take 25 lowest weights, then anything below threshold
    resampledPoints = []
    for point in samplePoints:
        #print(point.weight)
        if len(resampledPoints) < 25:
            resampledPoints.append(point)
        else:
            maxWeight = 0
            thisPoint = 0
            for addedPoint in resampledPoints:
                if addedPoint.weight > maxWeight:
                    maxWeight = addedPoint.weight
                    thisPoint = addedPoint
            if point.weight < maxWeight:
                resampledPoints.remove(thisPoint)
                resampledPoints.append(point)
            # elif point.weight < 0.1:
            #     resampledPoints.append(point)

    # Generate 400 new points grouped around our resampled particles
    finalParticles = []
    for point in resampledPoints:
        finalParticles.append(point)
        for i in range(4):
            newPoint = point.gaussianPoint(0.2)
            pointBody = Body(r=[AU * newPoint.position[0], AU * newPoint.position[1]], v=[0, 0], mass=0, measurement_noise=0)
            pointMeasure = pointBody.compute_gravity_magnitude(planets=solar_system.planets)
            newPoint.weight = abs(pointMeasure - objectMeasure)
            newPoint.measure = pointMeasure
            finalParticles.append(newPoint)

    # Step 4 - Fuzz particles
    # Generate some noise around particle

    # normal distribution around each particle???


    # Step 5 - Move based on satellite motion
    # ejected into circular counter-clockwide orbit around sun at position 0
    for point in finalParticles:
        angle = atan2(point.position[1], point.position[0])
        heading = angle * pi / 2
        point.direction = heading
        #point.move() # rotates counter-clockwise based on radius


    newTimer = Timer(timer.time+1, finalParticles)
    # Estimate
    # want to pick best guess
    minWeight = 1000
    minParticle = 0
    for particles in finalParticles:
        if particles.weight < minWeight:
            minWeight = particles.weight
            minParticle = particles

    ### Loop back to step 2 (assign important weights to new particles


    xy_estimate = (minParticle.position[0]*AU, minParticle.position[1]*AU)
    print(xy_estimate)

    # You may optionally also return a list of (x,y,h) points that you would like
    # the PLOT_PARTICLES=True visualizer to plot for visualization purposes.
    # If you include an optional third value, it will be plotted as the heading
    # of your particle.

    optional_points_to_plot = []
    for point in finalParticles:
        optional_points_to_plot.append((point.position[0]*AU, point.position[1]*AU, point.direction))
    #optional_points_to_plot = [(1*AU, 1*AU), (2*AU, 2*AU), (3*AU, 3*AU)]  # Sample (x,y) to plot
    #optional_points_to_plot = [(1*AU, 1*AU, 0.5), (2*AU, 2*AU, 1.8), (3*AU, 3*AU, 3.2)]  # (x,y,heading)

    return xy_estimate, newTimer, optional_points_to_plot


def next_angle(solar_system, gravimeter_measurement, other=None):
    """
    Gets the next angle at which to send out an sos message to the home planet,
    the last planet in the solar system's planet list.
    This is the function you will have to write for part B.

    The input parameters are exactly the same as for part A.

    :return:
        bearing: float. The absolute angle from the satellite to send an sos message
        estimate: Tuple[float, float]. The (x,y) estimate of the target satellite at the next timestep
        other: any. Any additional information you'd like to pass between invocations of this function
        optional_points_to_plot: List[Tuple[float, float, float]].
            A list of tuples like (x,y,h) to plot for the visualization
    """
    # At what angle to send an SOS message this timestep
    bearing = 0.0
    estimate = (110172640485.32968, -66967324464.19617)

    # You may optionally also return a list of (x,y) or (x,y,h) points that
    # you would like the PLOT_PARTICLES=True visualizer to plot.
    #
    # optional_points_to_plot = [ (1*AU,1*AU), (2*AU,2*AU), (3*AU,3*AU) ]  # Sample plot points
    # return bearing, estimate, other, optional_points_to_plot

    return bearing, estimate, other


def who_am_i():
    # Please specify your GT login ID in the whoami variable (ex: jsmith322).
    whoami = 'jgrc3'
    return whoami


# we will use the particle class
class Particle:
    mass_sun = 1.98847e30
    G = 6.6743e-11
    def __init__(self, particlePosition):
        self.position = particlePosition
        self.weight = 1000
        self.angle = atan2(particlePosition[1], particlePosition[0])
        self.measure = 0
        self.direction = 0

    # move particle --- counterclockwise motion around radius
    def move(self):
        radius = sqrt(self.position[0]**2 + self.position[1]**2)
        angle = atan2(self.position[1],self.position[0])

        vMag = sqrt(G * 1.98847e30)
        heading = angle * pi / 2
        self.direction = heading
        velocityX = vMag * cos(heading) / AU
        velocityY = vMag * sin(heading) / AU
        self.position = (self.position[0]+velocityX, self.position[1]+velocityY)
        # # velocity juts given in m/s --- need to calculate along radial path
        # arcMeasure = velocity / radius
        # newAngle = self.angle + arcMeasure
        # newPositionX = radius * math.cos(newAngle)
        # newPositionY = radius * math.sin(newAngle)
        # self.direction = (newPositionY - self.position[1]) / (newPositionX - self.position[0])
        # self.position = (newPositionX, newPositionY)

    def gaussianPoint(self, distance):
        newPositionX = self.position[0] + random.uniform(-1 * distance, distance)
        newPositionY = self.position[1] + random.uniform(-1 * distance, distance)
        return Particle((newPositionX, newPositionY))

# iteration
class Timer:
    def __init__(self, initialTime, samplePoints):
        self.time = initialTime
        self.particles = samplePoints