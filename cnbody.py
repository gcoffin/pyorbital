# requires numpy, tested on python 2.7
from __future__ import print_function

import cmath
import csv
import math
import random

import numpy as np

'''this module provides the classes to simulate a set of bodies subject to gravity.
All calculations are made in the plan, and positions are coded with complex number
to enable simple numpy calculations.

http://www.artcompsci.org/kali/vol/two_body_problem_2/ch11.html
'''

G = 6.6742E-11

COLL_TIME_MULTIPLIER = 0.01
VECTORTYPE = np.complex128


def vector(point):
    """helper to instantiate a point. In this implementation, converts
    a tuple to a complex"""
    if isinstance(point, complex):
        return point
    x, y = point
    return VECTORTYPE(x + 1j * y)


def inner(a, b):
    """inner product ax*b.x+a.y*b.y"""
    return a * np.conj(b)


def fformat(f):
    a = abs(f)
    if a > 10000 or a < 0.1:
        return '%.2g' % f
    if a > 100:
        return '%.0f' % f
    return '%.2f' % f


def cformat(c):
    return '[%s,%s]' % (fformat(c.real), fformat(c.imag))


class NBodySystem(object):
    """base implementation of N body problem with leadfrog integration
    in the plan (only x,y). Points are complex values."""

    def __init__(self, r, v, m, coll_time_multiplier):
        '''r: vector of position of the bodies at t=0, complex in m
        v: velocity, complex in (m/s)
        m: mass (kg)
        coll_time_multiplier: multiplier to the time step, typically 0.01
        All vectors of the same sizes.
        '''
        self.r = r
        self.v = v
        self.m = m
        self.N = len(m)
        self.opt_id = np.identity(self.N)
        self.opt_gmm = -G * np.multiply.outer(self.m, self.m)
        self.post_set()
        self.a, self.jk = self.acc_and_jerk(self.r, self.v)
        self.coll_time_multiplier = coll_time_multiplier

    # not tested yet.
    def add_body(self, r, v, m):
        self.r = np.append(self.r, r)
        self.v = np.append(self.v, v)
        self.m = np.append(self.m, m)
        self.N = len(m)
        self.opt_id = np.identity(self.N)
        self.opt_gmm = -G * np.multiply.outer(self.m, self.m)
        self.post_set()
        self.a, self.jk = self.acc_and_jerk(self.r, self.v)

    # hook to add post calculation. not used anymore
    def post_set(self):
        pass

    def acc_and_jerk(self, r, v):
        '''calculates the acceleration and jerk from positions and
        velocity of the bodies. Jerk is da/dt.'''
        rji = np.add.outer(r, -r)  # matrix of the relative positions
        vji = np.add.outer(v, -v)  # matrix of the relative velocitys
        r1 = np.absolute(rji)  # matrix of distances
        r2 = r1 * r1
        r3 = r1 * r2
        am = rji / r3 * self.m
        np.putmask(am, self.opt_id, 0.)
        a = -G * np.add.reduce(am, 1)
        jkm = -G * self.m / r3 * (vji - 3 * rji * inner(rji, vji) / r2)
        np.putmask(jkm, self.opt_id, 0.)
        jk = np.add.reduce(jkm, 1)
        return a, jk

    def step(self, t, dt):
        '''applies one step: calculates acceleration and jerk, then
        applies them to position and velocity'''
        old_r = self.r
        old_v = self.v
        old_a = self.a
        old_jk = self.jk
        # predictor
        dt2 = dt * dt
        r = self.r + self.v * dt + 0.5 * self.a * dt2 + self.jk * dt2 * dt / 6.
        v = self.v + self.a * dt + 0.5 * self.jk * dt2
        # applies with 2nd order
        self.a, self.jk = self.acc_and_jerk(r, v)
        self.v = old_v + (old_a + self.a) * (dt * .5) + (old_jk - self.jk) * dt2 / 12
        self.r = old_r + (old_v + self.v) * (dt * .5) + (old_a - self.a) * dt2 / 12
        self.post_set()

    def forward_fixed(self, starttime, endtime, dt=-1, steps=-1):
        '''moves from starttime to endtime.'''
        dT = float(endtime - starttime)
        dt = min(dt, dT)
        if steps == -1:
            truedt = dT / int(dT / dt)
            steps = int(dT / truedt + 0.5)
        else:
            truedt = dT / float(steps)
        for i in range(steps):
            self.step(starttime, truedt)
            starttime += truedt
        return steps, starttime

    def collision_time(self):
        '''Returns the min of (ri/vi), which is used to estimate the next dt'''
        # calculates the collision time -> review next dt
        rji = np.add.outer(self.r, -self.r)  # matrix of the relative positions
        vji = np.add.outer(self.v, -self.v)  # matrix of the relative velocitys
        v1 = np.absolute(vji)
        sel = np.nonzero(v1)
        return np.min(np.nan_to_num(np.absolute(rji)[sel] / v1[sel]))

    def energy(self):
        '''calculates the total energy (potential and kinetic) of the
        system'''
        rji = np.add.outer(self.r, -self.r)  # matrix of the relative positions
        r = np.absolute(rji)
        epotm = self.opt_gmm / r
        np.putmask(epotm, self.opt_id, 0.)
        epot = -0.5 * np.sum(epotm)
        ekin = 0.5 * np.add.reduce(self.m * inner(self.v, self.v))
        return (ekin + epot).real

    def barycenter(self):
        return self.r * self.m / np.sum(self.m)

    def speed_qty(self):
        return np.sum(self.m * self.v)

    def set_zero_speed_qty(self):
        self.v -= self.speed_qty() / np.sum(self.m)  # initialize total speed qty


class NBodyCommandSystem(NBodySystem):
    '''a n body system where we can apply a command'''

    def __init__(self, r, v, m, radius, coll_time_multiplier):
        self.acc_command = np.zeros(len(r), VECTORTYPE)
        self.radius = radius
        NBodySystem.__init__(self, r, v, m, coll_time_multiplier)
        self.sumradius = np.add.outer(self.radius, self.radius)
        np.putmask(self.sumradius, self.opt_id, 0)  # otherwise will collide planets with themselves

    def acc_and_jerk(self, r, v):
        a, jk = NBodySystem.acc_and_jerk(self, r, v)
        a += self.acc_command
        return a, jk

    def post_step(self, dt):
        pass

    def step(self, t, dt):
        NBodySystem.step(self, t, dt)
        self.post_step(dt)

    def check_collide(self):
        '''checks if the distance between 2 bodies is smaller than the
        sum of their radiuses. This does not take into account
        motion, but with the collision_time this should not be a
        problem'''

        # self.sumradius=abs(np.add.outer(self.radius,self.radius))
        # np.putmask(self.sumradius,self.opt_id,0) #otherwise will collide planets with themselves
        rji = np.add.outer(self.r, -self.r)  # matrix of the relative positions
        distmatrix = abs(rji) - self.sumradius
        result = []
        close_bodies = distmatrix < 0
        if np.sometrue(close_bodies):
            vji = np.add.outer(self.v, -self.v)  # matrix of the relative velocitys
            # they are close and going towards each other
            towards = -inner(vji, rji)
            collision = np.logical_and(close_bodies, towards > 0)
            if np.any(collision):
                return [(i, j, towards[i, j]) for i, j in zip(*np.nonzero(collision)) if i < j]
        return []

    def closest(self, idx, ignorelist=()):
        '''gets the closest body from body (indexes in the list)'''
        ignore = set(ignorelist)
        order = (np.absolute(self.r - np.ones(self.r.shape) * self.r[idx]) - self.radius).argsort()
        for idx in order[1:]:  # first is the body itself
            if idx not in ignore:
                return idx
        return None

    def attractor(self, idx):
        r = self.r - self.r[idx]
        a = np.absolute(self.m / (r * r))
        a[idx] = 0
        return a.argmax()


class BoundSystemDescriptor(object):
    def __init__(self, name):
        self.name = name

    def __get__(self, instance, owner=None):
        return getattr(instance._system, self.name)[instance._index]

    def __set__(self, instance, value):
        getattr(instance._system, self.name)[instance._index] = value


class Body(object):
    '''An element subject to gravitational force, this is a proxy to
    the elements stored in the system. First (1) instantiate it, then
    (2) use all bodies to create a system, then (3) bind all bodies to
    the system with the method bind. '''

    def __init__(self, position, velocity, mass):
        self.initialposition = vector(position)
        self.initialvelocity = vector(velocity)
        self.mass = mass
        self._system = None
        self._index = -1

    def info(self):
        return "<body (%.3e,%.3e) (%.3e,%.3e) %.3e>" % (self.r.real, self.r.imag, self.v.real, self.v.imag, self.mass)

    def bind(self, system, index):
        '''binds this body to a system'''
        self._system = system
        self._index = index

    @property
    def idx(self):
        return self._index

    v = BoundSystemDescriptor('v')
    a = BoundSystemDescriptor('a')
    r = BoundSystemDescriptor('r')
    m = BoundSystemDescriptor('m')
    jk = BoundSystemDescriptor('jk')


class SpaceBody(Body):
    '''a body with a radius and an atmosphere for use in NBodyCommandSystem'''

    def __init__(self, name, position, velocity, mass, radius, parent, atm, **attrs):
        Body.__init__(self, position, velocity, mass)
        self.initialradius = radius
        self.name = name
        self.parent = parent
        self.atm = atm
        for i, j in attrs.items():
            setattr(self, i, j)

    def __str__(self):
        return self.name

    def __repr__(self):
        return "<SpaceBody %s>" % self.name

    @property
    def radius(self):
        return self._system.radius[self._index]

    @radius.setter
    def radius(self, value):
        self._system.radius[self._index] = value
        self._system.sumradius = np.add.outer(self._system.radius, self._system.radius)

    @property
    def closest(self):
        return self._system.closest_body(self)

    @property
    def attractor(self):
        return self._system.bodylist[self._system.attractor(self.idx)]

    @property
    def distance_to_closest(self):
        return abs(self.r - self.closest.r)

    @property
    def rel_coord(self):
        class Coord:
            def __init__(self, body):
                self.body = body
                self.parent = self.body.closest
                r = body.r - body.closest.r
                self.u = (r / abs(r)).conjugate()

            def __getattr__(self, attr):
                return (getattr(self.body, attr) - getattr(self.parent, attr)) * self.u

        return Coord(self)

    # the following calculation are relative to another body assuming a 2-body system
    def barycenter(self, body):
        return (body.m * body.r + self.m * self.r) / (body.m + self.m)

    def barycenter_vel(self, body):
        return (body.m * body.v + self.m * self.v) / (body.m + self.m)

    def rel_energy(self, body):
        r = self.r - body.r
        v = self.v - body.v
        mm = self.m * body.m
        epot = -G * mm / np.absolute(r)
        ekin = 0.5 * self.m * np.absolute(v) ** 2
        return (epot + ekin)

    def escape_speed(self, body):
        return math.sqrt(2 * G * body.m / np.absolute(body.r - self.r))

    def eccentricity(self, body):
        rel_pos = self.r - self.barycenter(body)
        rel_velocity = self.v - self.barycenter_vel(body)
        red_mass = (self.m + body.m) / (self.m * body.m)
        ang_momentum = (rel_velocity.conjugate() * rel_pos).imag  # cross vector
        ecc2 = max(0, 1 + 2 * self.rel_energy(body) * ang_momentum ** 2 * red_mass / ((body.m * G) ** 2))
        return math.sqrt(ecc2)


class SpaceShip(SpaceBody):
    '''a body with acceleration command for use in NBodyCommandSystem'''

    acc_command = BoundSystemDescriptor('acc_command')

    @property
    def acc_command_polar(self):
        '''get the command in polar coordinate: r,phi in degres'''
        a, phi = cmath.polar(self.acc_command)
        return a, phi / math.pi * 180.

    @acc_command_polar.setter
    def acc_command_polar(self, command):
        (acc, angle) = command
        self.acc_command = cmath.rect(acc, angle / 180. * math.pi)
        f = self.rel_coord


class BodyListSystem(NBodySystem):
    def __init__(self, bodylist, coll_time_multiplier):
        self.bodylist = bodylist
        m = np.array([i.mass for i in bodylist], np.float64)
        r = np.array([i.initialposition for i in bodylist], VECTORTYPE)
        v = np.array([i.initialvelocity for i in bodylist], VECTORTYPE)
        radius = np.array([i.initialradius for i in bodylist], np.float64)
        for i, b in enumerate(self.bodylist):
            b.bind(self, i)
        NBodySystem.__init__(self, r, v, m, coll_time_multiplier)

    @property
    def parents(self):
        return [body for body in self.bodylist if body.parent is None]

    def closest_body(self, body, *ignore):
        idx = self.closest(body.idx, [b.idx for b in ignore])
        return self.bodylist[idx]

    def check_body_collide(self):
        return [(self.bodylist[i], self.bodylist[j], t) for i, j, t in self.check_collide() if i < j]


class BodyListCommandSystem(NBodyCommandSystem, BodyListSystem):
    def __init__(self, bodylist, coll_time_multiplier):
        self.bodylist = bodylist
        m = np.array([i.mass for i in bodylist], np.float64)
        r = np.array([i.initialposition for i in bodylist], VECTORTYPE)
        v = np.array([i.initialvelocity for i in bodylist], VECTORTYPE)
        radius = np.array([i.initialradius for i in bodylist], np.float64)
        for i, b in enumerate(self.bodylist):
            b.bind(self, i)
        NBodyCommandSystem.__init__(self, r, v, m, radius, coll_time_multiplier)


class rnd:
    def __init__(self):
        self.vals = {}
        self.last = random.random()

    def __getitem__(self, i):
        if not i in self.vals:
            self.value[i] = random.random()
        return self.value[i]

    def rnd(self):
        self.last = random.random()
        return self.last


def get_context():
    '''functions for use in the formulas in the configuration file'''
    result = {}
    import math, cmath, random
    result.update(math.__dict__)
    result.update(cmath.__dict__)
    result.update(random.__dict__)
    result['rnd'] = rnd()
    return result


formula_context = get_context()


class NamedSystemBuilder(object):
    '''object to build a system from a csv file'''

    def __init__(self, *planets):
        self.bodylist = []
        self.bodymap = {}
        self.add(*planets)

    def getsystem(self, coll_time_multiplier):
        return BodyListSystem(self.bodylist, coll_time_multiplier)

    def makeplanet(self, name, distance, phaseindeg, direction, radius, mass, parentname, atm=False, cls=SpaceBody,
                   **attrs):
        '''Creates a planet given the parameters. The parameters may be strings or formulas to eval.
        name: name of the body. Doesn't have to be unique
        distance: distance to the parent in m (the star for a planet, a placet for a satellite, etc) or (0,0)
        phaseindeg: phase in degree.
        direction: velocity of the body compared to the velocity required to stay in orbit. 
           Enter 1 for a circle orbit, 1j to go towards the parent
        radius: radius in m
        mass: mass in kg
        parentname: name of the parent
        atm: true is there is an atmosphere
        cls: class to implement the body. SpaceBody by default, Spaceship for the spaceship
        attrs: additional parameters directly passed to the constructor of cls.
        '''
        phase = phaseindeg / 180. * math.pi
        if parentname == '--' or parentname == name:
            parent = None
        else:
            parent = self.bodymap[parentname]

        if parent:
            distance += parent.initialradius
        # find the velocity according to parameter "direction"
        u = cmath.exp(1j * phase)
        if parent is None:
            velocity = 0j
            position = distance * u
        else:
            position = distance * u + parent.initialposition
            gravityfromparent = G * parent.mass / distance ** 2
            velocity = math.sqrt(distance * gravityfromparent) * u * 1j * direction + parent.initialvelocity
        attr = {}
        if atm:
            for i, j in attrs.items():
                attr[i] = float(j)

        return cls(name, position, velocity, mass, radius, parent, atm, **attr)

    def makeplanet_eval(self, name, distance, phaseindeg, direction, radius, mass, parentname, atm='False',
                        cls=SpaceBody, **attrs):
        global formula_context
        name = name.strip()
        parentname = parentname.strip()
        distance = eval(distance, formula_context)
        phaseindeg = eval(phaseindeg, formula_context)
        direction = eval(direction, formula_context)
        radius = eval(radius, formula_context)
        mass = eval(mass, formula_context)
        atm = atm.strip().lower() == 'true'
        return self.makeplanet(name, distance, phaseindeg, direction, radius, mass, parentname, atm, cls, **attrs)

    def readfile(self, filename):
        '''read csv file'''
        reader = csv.DictReader(open(filename), skipinitialspace=True)
        for j in reader:
            # print (j)
            if j['name'].startswith("#"):
                continue
            if j['name'].strip() == 'Spaceship':
                j['cls'] = SpaceShip
            planet = self.makeplanet_eval(**j)
            self.add(planet)

    def execfile(self, filename):
        context = formula_context.copy()

        def make(*args, **kwargs):
            self.add(self.makeplanet(*args, **kwargs))

        context['make'] = make
        context['SpaceShip'] = SpaceShip
        exec(open(filename).read(), context)

    def add(self, *planets):
        for planet in planets:
            self.bodylist.append(planet)
            self.bodymap[planet.name] = planet

    def getbodylist(self):
        return self.bodylist


def main(filename):
    b = NamedSystemBuilder()
    b.readfile(filename)
    s = b.getsystem(0.1)

    rows = (('time', 't'),
            ('name', 'i.name'),
            ('x', 'i.r.real/1e6'),
            ('y', 'i.r.imag/1e6'),
            ('ditance', 'abs(i.r-getparentposition(i.parent))'),
            ('velocity', 'abs(i.v-getparentposition(i.parent))')
            )

    t = 0
    TOTAL_T = 100000
    N = 10
    dt = TOTAL_T / N
    en0 = s.energy()
    for i in range(N):
        steps, t = s.forward_fixed(t, t + dt, dt=1000)
    en1 = s.energy()
    print("Error on energy", (en1 - en0) / en0)


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        if sys.argv[1] == '-p':
            print('profiling')
            import cProfile

            cProfile.run('main("planet.csv")', 'profile')
        elif sys.argv[1][-4:] == '.csv':
            main(sys.argv[1])
    else:
        main('missionmoon.csv')
