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
    timer = other
    objectMeasure = gravimeter_measurement
    samplePoints = []
    # for i in range(1000):
    #     randNumX = random.random() * 8 - 4
    #     randNumY = random.random() * 8 - 4
    #     thisPoint = Particle((randNumX, randNumY))
    #     samplePoints.append(thisPoint)
    if not timer:
        # want to initialize 1000 random points near object
        for i in range(10000):
            randNumX = random.random() * 8 - 4
            randNumY = random.random() * 8 - 4
            thisPoint = Particle((randNumX, randNumY))
            samplePoints.append(thisPoint)
        # avgWeight = 0
        # for i in range(500):
        #     randNumX = random.random() * 8 - 4
        #     randNumY = random.random() * 8 - 4
        #     pointBody = Body(r=[AU * randNumX, AU * randNumY], v=[0, 0], mass=0, measurement_noise=00)
        #     pointMeasure = pointBody.compute_gravity_magnitude(planets=solar_system.planets)
        #     avgWeight += 1 / abs(pointMeasure - objectMeasure)  # take reciprocal s.t. weight large
        # avgWeight /= 500
        # for i in range(5000):
        #     weight = 0
        #     x = 0
        #     y = 0
        #     while weight < avgWeight * 2:
        #         randNumX = random.random() * 8 - 4
        #         randNumY = random.random() * 8 - 4
        #         pointBody = Body(r=[AU * randNumX, AU * randNumY], v=[0, 0], mass=0, measurement_noise=00)
        #         pointMeasure = pointBody.compute_gravity_magnitude(planets=solar_system.planets)
        #         weight = 1 / abs(pointMeasure - objectMeasure)  # take reciprocal s.t. weight large
        #     samplePoints.append(Particle((randNumX, randNumY)))
        time = 0
    else: # initialize as past set
        for i in range(len(timer.particles)):
            samplePoints.append(Particle(timer.particles[i].position))
        time = timer.time + 1

    # # Step 2 - Assign importance weights - how likely particle is to be object - assign based on gravimeter
    # # first need gravimeter measure
    # need to assign score to each point -- closer to 0, the closer the score
    for point in samplePoints:
        pointBody = Body(r=[AU * point.position[0], AU * point.position[1]], v=[0, 0], mass=0, measurement_noise=00)
        pointMeasure = pointBody.compute_gravity_magnitude(planets=solar_system.planets)
        ### weight needs to be gaussian!!!
        #point.weight = 1 / abs(pointMeasure - objectMeasure)
        point.weight = Gaussian(pointMeasure, 0.1, objectMeasure)  # take reciprocal s.t. weight large
        #point.measure = pointMeasure
        #print(pointMeasure)
    #
    #
    #
    # Step 3 - Resample - generate N new points based on importance weights of last set
    # keep the 25 best particles??? -- take 25 lowest weights, then anything below threshold
    # We need to use importance wheel
    # first get normalized weight
    totalWeight = 0
    for point in samplePoints:
        totalWeight += point.weight

    # assign norm weights
    norms = []
    for point in samplePoints:
        weight = point.weight / totalWeight
        point.normWeight = weight
        norms.append(point.normWeight)

    #print(norms)
    # resampling wheel
    N = 100
    resampledPoints = []
    index = int(random.random() * N)
    beta = 0.0
    mw = max(norms)
    for i in range(N):
        beta += random.random() * 2.0 * mw
        while norms[index] < beta:
            beta = beta - norms[index]
            index = (index + 1) % N
        resampledPoints.append(samplePoints[index])



    # Step 4 - Fuzz particles
    # Generate some noise around particle
    # normal distribution around each particle???
    # Generate 400 new points grouped around our resampled particles
    finalParticles = []
    # order points by weight
    resampledPoints.sort(key=lambda x: x.normWeight, reverse=True)

    for i in range(len(resampledPoints)):
        numberFuzz = int(10 / (int(i/10 + 1)))
        finalParticles.append(resampledPoints[i])
        for j in range(numberFuzz):
            # newPoint = resampledPoints[i].gaussianPoint(0.2)
            # pointBody = Body(r=[AU * newPoint.position[0], AU * newPoint.position[1]], v=[0, 0], mass=0,
            #                  measurement_noise=0)
            # pointMeasure = pointBody.compute_gravity_magnitude(planets=solar_system.planets)
            # newPoint.weight = 1 / abs(pointMeasure - objectMeasure)
            # newPoint.measure = pointMeasure
            # finalParticles.append(newPoint)

            # add radius fuzz
            radius = sqrt(resampledPoints[i].position[0]**2 + resampledPoints[i].position[1]**2)
            theta = atan2(resampledPoints[i].position[1], resampledPoints[i].position[0]) + random.uniform(-1*pi/32, pi/32)
            x = radius * cos(theta)
            y = radius * sin(theta)
            radiusRand = Particle((x,y))
            newPoint = radiusRand.gaussianPoint(0.1)
            radBody = Body(r=[AU * newPoint.position[0], AU * newPoint.position[1]], v=[0, 0], mass=0, measurement_noise=0)
            radMeasure = radBody.compute_gravity_magnitude(planets=solar_system.planets)
            #newPoint.weight = 1 / abs(radMeasure - objectMeasure)
            newPoint.weight = Gaussian(radMeasure, 0.1, objectMeasure)
            finalParticles.append(newPoint)


    # # Step 5 - Move based on satellite motion
    # # ejected into circular counter-clockwide orbit around sun at position 0
    movedPoints = []
    for point in finalParticles:
        newPoint = point.moveBody() # rotates counter-clockwise based on radius
        movedPoints.append(newPoint)

    movedPoints.sort(key=lambda x: x.weight, reverse=True)

    maxParticle = movedPoints[0]
    # maxX = []
    # maxY = []
    # indexWeights = []
    # sumWeight = 0
    # for point in movedPoints[:10]:
    #     indexWeights.append(point.weight)
    #     maxX.append(point.position[0])
    #     maxY.append(point.position[1])
    # for i in range(len(indexWeights)):
    #     sumWeight += indexWeights[i]
    # for i in range(len(indexWeights)):
    #     indexWeights[i] /= sumWeight

    # posX = 0
    # posY = 0
    # for i in range(len(indexWeights)):
    #     posX += maxX[i] * indexWeights[i]
    #     posY += maxY[i] * indexWeights[i]

    # ### Loop back to step 2 (assign important weights to new particles

    xy_estimate = (maxParticle.position[0]*AU, maxParticle.position[1]*AU)#, movedPart.direction)

    # #print(xy_estimate)
    #
    # # You may optionally also return a list of (x,y,h) points that you would like
    # # the PLOT_PARTICLES=True visualizer to plot for visualization purposes.
    # # If you include an optional third value, it will be plotted as the heading
    # # of your particle.
    #
    optional_points_to_plot = []
    for point in movedPoints[:100]:
        optional_points_to_plot.append((point.position[0]*AU, point.position[1]*AU, point.direction))

    timer = Timer(time, movedPoints)

    return xy_estimate, timer, optional_points_to_plot


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
    #estimate = (110172640485.32968, -66967324464.19617)

    # You may optionally also return a list of (x,y) or (x,y,h) points that
    # you would like the PLOT_PARTICLES=True visualizer to plot.
    #
    # optional_points_to_plot = [ (1*AU,1*AU), (2*AU,2*AU), (3*AU,3*AU) ]  # Sample plot points
    # return bearing, estimate, other, optional_points_to_plot
    position, other, optional_points_to_plot = estimate_next_pos(solar_system, gravimeter_measurement, other)

    estimate = position
    lastPlanet = solar_system.planets[-1]

    # need to move last planet one ahead
    lastPlanetMove = solar_system.move_body(lastPlanet)
    lastPlanetPos = (lastPlanetMove.r[0], lastPlanetMove.r[1])

    # now get bearing from estimate to lastPlanetPos
    bearing = atan2(lastPlanetPos[1]-estimate[1], lastPlanetPos[0]-estimate[0])


    return bearing, estimate, other

def Gaussian(mu, sigma, x):
    # calculates the probability of x for 1-dim Gaussian with mean mu and var. sigma
    return exp(- ((mu - x) ** 2) / (sigma ** 2) / 2.0) / sqrt(2.0 * pi * (sigma ** 2))


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
        self.direction = atan2(particlePosition[1], particlePosition[0]) + pi/2
        self.normWeight = 0

    # move particle --- counterclockwise motion around radius
    def move(self):
        #pointBody = Body(r=[AU * self.position[0], AU * self.position[1]], v=[0, 0], mass=0, measurement_noise=0)
        #movedBody = pointBody.moveBody(planets)
        # movedPart = Particle((movedBody.r[0], movedBody.r[1]))
        # movedPart.weight = self.weight
        radius = sqrt(self.position[0]**2 + self.position[1]**2)
        heading = self.direction
        vMag = sqrt(G * 1.98847e30)
        # velocity_magnitude = vMag / sqrt(radius * AU)  # m/s
        # heading = heading + pi / 2  # perpendicular
        # velocity_x = velocity_magnitude * cos(heading) / AU
        # velocity_y = velocity_magnitude * sin(heading) / AU
        #vMag2 = sqrt(G * 1.98847e30 / (radius * AU))
        arcMeasure = vMag / (radius * AU)
        newAngle = self.angle + arcMeasure
        newPositionX = radius * math.cos(newAngle)
        newPositionY = radius * math.sin(newAngle)
        #heading = atan2(self.position[1], self.position[0]) * pi/2
        #velocityX = vMag * math.cos(heading) / radius
        #velocityY = vMag * math.cos(heading) / radius
        #print(heading, velocityX, velocityY)
        #position = (self.position[0]+velocity_x, self.position[1]+velocity_y)
        position = (newPositionX, newPositionY)
        movedPart = Particle(position)
        movedPart.weight = self.weight
        return movedPart

    def moveBody(self):
        x = self.position[0]
        y = self.position[1]
        radius_body = sqrt(x ** 2 + y ** 2)
        angle = atan2(y, x)
        velocity_magnitude = sqrt(G * 1.98847e30 / radius_body)  # m/s
        heading = angle + pi / 2  # perpendicular
        velocity_x = velocity_magnitude * cos(heading)
        velocity_y = velocity_magnitude * sin(heading)
        newX = x + (velocity_x / AU)
        newY = y + (velocity_y / AU)
        movedParticle = Particle((newX, newY))
        movedParticle.weight = self.weight
        return movedParticle

    def gaussianPoint(self, distance):
        newPositionX = self.position[0] + random.uniform(-1 * distance, distance)
        newPositionY = self.position[1] + random.uniform(-1 * distance, distance)
        return Particle((newPositionX, newPositionY))


# iteration
class Timer:
    def __init__(self, initialTime, samplePoints):
        self.time = initialTime
        self.particles = samplePoints
