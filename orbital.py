# requires: numpy and pygame

import cmath
import math
import os
import re
import time

import numpy
import pygame
from pygame.locals import *

import cnbody
import gameloop as gl
import loader

ACC_MULT0 = 8  # initial spaceship thrust - 2g
FUEL = 50000.  # initial volume of fuel = DeltaV

ZOOM_QUANTUM = 2
COLL_TIME_MULTIPLIER = .01  # calculates desired calculation time step according to typical system time
ZOOM_EXP = 0.55  # so that a constant acceleration produces the same visual effect at all scales
ZOOM_TIME_MULT = 0.1
ROTATION_DEG_BY_SEC = 180
IMPACT = 1e40
COEFFRESTITUTION = 0.9

AUTOZOOM_TIMESTEP_MULT = 8
AUTOZOOM_QUANTIZE = 3
KEYPRESSED_EVENT = 42
TIME_STEP = 60

MAX_IMAGE_SIZE = 4000.  # maximum image size allowed by pygame

FONT_NAME = None  # 'vera.ttf' #'digitalism.ttf' #


def inner(u, v):
    """inner product"""
    return (u * v.conjugate()).real


class CrashException(Exception):
    pass


class GameNBodySystem(cnbody.BodyListCommandSystem):
    def __init__(self, bodylist, coll_time_multiplier):
        self.spaceship_body = None
        for i in bodylist:
            if i.name == 'Spaceship':
                self.spaceship_body = i
        cnbody.BodyListCommandSystem.__init__(self, bodylist, coll_time_multiplier)

    def acc_and_jerk(self, r, v):
        a, jk = cnbody.BodyListCommandSystem.acc_and_jerk(self, r, v)
        # atmospherical
        if not self.spaceship_body:
            return a, jk
        si = self.spaceship_body._index
        ci = self.closest(si)
        r = abs(self.r[si] - self.r[ci])
        if r < self.radius[ci] * 1.5:
            x = self.bodylist[ci]
            if x.atm:
                h = max(0, r - self.radius[ci])
                p = x.surf_p * math.exp(-x.surf_g * x.molar_mass * h / (8.31432 * x.surf_temp))
                if p > 1e3:
                    v = (self.spaceship_body.v - x.v)
                    friction = -v * abs(v) * p * 5e-10
                    if friction > 200:  # 20g: boom
                        raise CrashException('')
                    a[self.spaceship_body._index] += friction
        return a, jk

    def post_step(self, dt):
        collision = self.check_body_collide()
        if collision:
            for i, j, t in collision:
                if self.spaceship_body in (i, j):
                    i, j = (i, j) if i == self.spaceship_body else (j, i)
                    if abs(i.v - j.v) > 300:  # > 1000 km/h, boom!
                        raise CrashException("crashed %f" % (abs(i.v - j.v)))
                    if t > 0:  # landed
                        i.v = j.v
                        i.r = j.r + (i.r - j.r) / abs(i.r - j.r) * (i.radius + j.radius)
                else:
                    self.bounces(i, j)
            self.post_set()

    def bounces(self, body, body2):
        """2 bodies bounce on each other"""
        coeffrestitution = getattr(body, 'coeffrestitutison', COEFFRESTITUTION) * getattr(body2, 'coeffrestitution',
                                                                                          COEFFRESTITUTION)
        r = body.r - body2.r
        d = abs(r)
        u = r / d  # unit vector
        n = r * 1j  # normal vector
        tm = (body.mass + body2.mass)
        vb = (body.mass * body.v + body2.mass * body2.v) / tm  # velocity of the barycenter
        sradii = body.radius + body2.radius
        if d < sradii:
            # too close? adjust so they are in contact
            # consider strong force instead?
            rb = (body.mass * body.r + body2.mass * body2.r) / tm
            body.r = rb + u * sradii * body2.mass / tm
            body2.r = rb - u * sradii * body.mass / tm
        v1 = (body.v - vb) * coeffrestitution  # keep speed quantity
        v2 = (body2.v - vb) * coeffrestitution
        body.v = (v1 / n).conjugate() * n + vb  # bounces
        body2.v = (v2 / n).conjugate() * n + vb


class NamedPosition(object):
    def __init__(self, name, coords):
        self.r = cnbody.vector(coords)
        self.name = name

    def __repr__(self):
        return "<NamedPosition %s>" % self.name


class OnSurface(NamedPosition):
    def __init__(self, planet, point):
        self.planet = planet
        self.point = point
        self.name = planet.name

    @property
    def r(self):
        u = self.point.r - self.planet.r
        return self.planet.r + self.planet.radius * u / abs(u)

    def __repr__(self):
        return "<OnSurface projects %s on %s>" % (self.point, self.planet)


class Barycenter(NamedPosition):
    def __init__(self, system):
        super().__init__('',None)
        self.system = system
        self.m = numpy.sum(self.system.m)

    @property
    def r(self):
        return numpy.sum(self.system.m * self.system.r) / self.m


DEFAULT_POS = NamedPosition("--", (0, 0))


class ScreenTransform(object):
    """converts physical position to screen position"""

    def __init__(self, wh, xy, origin=None, viewradius=-1, scale=-1):
        """width,height:size of the surface
        x,y: screen position of origin
        """
        self.scale = None
        (width, height) = wh
        self.sx_position, self.sy_origin = xy
        self.set_frame_size((width, height))
        if scale > 0:
            self.set_scale(scale)
        else:
            self.set_viewradius(viewradius)
        if origin is None:
            origin = DEFAULT_POS
        self.origin = origin

    def set_frame_size(self, wh):
        (width, height) = wh
        self.width = width
        self.height = height
        self.screenradius = min(self.width, self.height) / 2

    def set_viewradius(self, viewradius):
        """set scale according to the to be viewed on the frame"""
        self.scale = self.screenradius / viewradius

    def set_scale(self, scale):
        self.scale = scale
        return self.scale

    def convert(self, pos):
        """converts a real position to pixel coordinates"""
        pos -= self.origin.r
        return (self.sx_position + int(pos.real * self.scale),
                self.height - self.sy_origin - int(pos.imag * self.scale))

    def backconvert(self, spos):
        """converts pixel coordinates to real position"""
        x, y = spos
        return (x - self.sx_position + 1j * (y - self.sy_position)) / self.scale + self.origin.r

    def screen_rect(self, rect):
        return self.backconvert(rect.topleft), self.backconvert(rect.bottomright)

    def set_pos(self, sprite):
        """positions a sprite on the screen"""
        x, y = self.convert(sprite.body.r)
        sprite.set_pos(x, y)

    def safe_set_pos(self, sprite):
        """positions a sprite on the screen, except if the sprite is
        outside the viewing area. Returns True if the sprite is
        visible"""
        x, y = self.convert(sprite.body.r)
        radius = int(sprite.body.radius * self.scale)
        if x + radius < 0 or y + radius < 0 or x - radius > self.width or y - radius > self.height:
            result = False
        else:
            # sprite.rect.topleft=x-sprite.rect.width/2,y-sprite.rect.height/2 #topleft center
            sprite.set_pos(x, y)
            result = True
        return result


class BodySprite(pygame.sprite.DirtySprite):
    """representation of a body on a frame"""

    def __init__(self, body, frame, imagename, imagesize=-1):
        pygame.sprite.DirtySprite.__init__(self)
        self.body = body
        self.frame = frame
        # Intialize Sprite, which is our base class
        self.baseimage = loader.load_image(imagename, True)
        self.baseimagerect = self.baseimage.get_rect()
        if imagesize == -1:
            imagesize = self.baseimagerect.width
        self.baseimage_actradius = imagesize
        # self.rect = self.image.get_rect()
        self.size = None
        self.angle = self.body_angle = 0
        # redo the math: scale so that the scaled image won't be above a certain size (1000px for instance)
        self.maxscale = float(MAX_IMAGE_SIZE * self.baseimage_actradius) / (2 * self.body.radius * imagesize)
        self.min_size = 6

    def set_image(self, scale, angle):
        apply_rotation = self.angle != angle
        size = max(self.min_size, 2 * self.body.radius * scale)
        if self.size != size:
            ratio = size / float(self.baseimage_actradius)
            self.image = self.base_rotimage = pygame.transform.scale(self.baseimage, (
            int(ratio * self.baseimagerect.width), int(ratio * self.baseimagerect.height)))
            self.size = size
            apply_rotation = True
        if apply_rotation:
            self.image = pygame.transform.rotate(self.base_rotimage, angle)
            self.angle = angle
            self.rect = self.image.get_rect()

    def update(self):
        self.visible = self.frame.safe_set_pos(self)
        self.dirty = True

    def rescale(self, scale):
        self.set_image(scale, self.body_angle)
        self.scale = scale

    def set_pos(self, x, y):
        self.rect.center = x, y


class SpaceShipSprite(BodySprite):
    def __init__(self, spaceship, frame, imageidle, imageforward, imagebackward, imageexplosion):
        BodySprite.__init__(self, spaceship, frame, imageidle)
        self.images = {-1: loader.load_image(imagebackward, True),
                       0: loader.load_image(imageidle, True),
                       1: loader.load_image(imageforward, True)}
        self.explosion = loader.load_image(imageexplosion)
        self.user_command_acc = 0.
        self.prev_body_angle = -1.
        self.body_angle = 0.
        self.prev_user_command_acc = 0
        self.rotate = True
        self.min_size = 24
        self.crashed = -1

    def crash(self):
        self.image = self.explosion
        if self.crashed == -1: self.crashed = 50

    def update(self):
        """check if the spaceship image needs change (ignite/stop fire, rotate)"""
        if self.crashed == -1:
            if self.user_command_acc < 0:
                i = -1
            elif self.user_command_acc > 0:
                i = 1
            else:
                i = 0
            if i != self.prev_user_command_acc or self.body_angle != self.prev_body_angle:
                self.baseimage = self.images[i]
                self.size = None
                body_angle = self.angle = 0
                if self.rotate: body_angle = self.body_angle
                self.set_image(self.scale, body_angle)
                self.prev_user_command_acc = self.user_command_acc
                self.prev_body_angle = body_angle
        else:
            self.crashed = max(0, self.crashed - 1)
        BodySprite.update(self)
        if self.crashed == 0:
            self.visible = False


def format_distance(num):
    if num > 1E9:
        return "%.2f M km" % (num / 1e9)
    elif num > 1E6:
        return "%d %03d km" % (int(num / 1e6), int(num / 1e3 % 1e3))
    else:
        return "%.0f km" % (num / 1e3)


def format_speed(speed):
    return "%.0f km/h" % (speed * 3.6)


class TextDisplay(pygame.sprite.Sprite):
    def __init__(self, rect, header,
                 backgroundcolor=(210, 210, 210),
                 fontcolor=(0, 0, 0),
                 adaptheight=True,
                 refreshfreq=10,
                 margin=3):
        pygame.sprite.Sprite.__init__(self)
        self.updatecount = 0
        self.refreshfreq = refreshfreq
        # initializes header
        self.header = []
        self.fonts = {}
        self.split = 0
        height = 3
        for name, fontsize in header:
            if fontsize not in self.fonts:
                self.fonts[fontsize] = loader.load_font(FONT_NAME, fontsize)
            w, h = self.fonts[fontsize].size(name)
            self.split = max(w, self.split)
            self.header.append((name, fontsize, height + margin))
            height += h + margin * 2
        if adaptheight:
            rect.height = height + margin
        self.split += 8
        self.left = rect.left
        # initialize background and colors
        self.fontcolor = fontcolor
        self.image = pygame.Surface(rect.size)
        self.backgroundcolor = backgroundcolor
        self.image.fill(self.backgroundcolor)
        self.rect = self.image.get_rect()
        self.rect.topleft = rect.topleft

    def update(self):
        self.updatecount -= 1
        if self.updatecount > 0:
            return
        self.updatecount = self.refreshfreq

        self.image.fill(self.backgroundcolor)
        self.render()

    def render(self):
        for (name, fontsize, y) in self.header:
            self.image.blit(self.fonts[fontsize].render(name, 1, self.fontcolor), (self.left, y))


class DashboardDisplay(TextDisplay):
    def __init__(self, rect, text,
                 backgroundcolor=(210, 210, 210),
                 fontcolor=(0, 0, 0),
                 adaptheight=True,
                 refreshfreq=20,
                 margin=3):
        TextDisplay.__init__(self, rect, text, backgroundcolor, fontcolor, adaptheight, refreshfreq, margin)
        self.values = {}

    def render(self):
        for (name, fontsize, y) in self.header:
            self.image.blit(self.fonts[fontsize].render(name, 1, self.fontcolor), (self.left, y))
            self.image.blit(self.fonts[fontsize].render(self.values.get(name, '??'), 1, self.fontcolor),
                            (self.left + self.split, y))


class SimpleFrame(object):
    """a frame is a surface and sprites"""

    def __init__(self, name, screensurface, backgroundimage, backgroundcolor=None, spritegroup=None):
        self.name = name
        if spritegroup is None:
            spritegroup = pygame.sprite.RenderUpdates()
        self.spritegroup = spritegroup
        self.screensurface = screensurface
        self.backgroundimage = pygame.transform.scale(backgroundimage, screensurface.get_rect().size)
        self.screensurface.blit(self.backgroundimage, (0, 0))

    def add(self, sprite):
        self.spritegroup.add(sprite)

    def clear(self):
        self.spritegroup.clear(self.screensurface, self.backgroundimage)

    def draw(self):
        self.spritegroup.draw(self.screensurface)

    def update(self):
        self.spritegroup.update()

    def onceinawhile(self):
        pass


class Frame(ScreenTransform):
    """a frame that transforms object position into screen position"""

    def __init__(self, name, screensurface, origin, viewradius, backgroundimage, backgroundcolor=None,
                 spritegroup=None):
        self.name = name
        self.prevscale = -1
        self.rect = screensurface.get_rect()
        if spritegroup is None:
            spritegroup = pygame.sprite.LayeredDirty()
        self.spritegroup = spritegroup

        ScreenTransform.__init__(self, self.rect.size, self.rect.center, viewradius=viewradius, origin=origin)

        self.screensurface = screensurface
        self.backgroundimage = pygame.transform.scale(backgroundimage, screensurface.get_rect().size)
        self.screensurface.blit(self.backgroundimage, (0, 0))

    def set_viewradius(self, viewradius):
        self.set_scale(self.screenradius / viewradius)

    def is_visible(self, body, scale):
        return abs(body.r - self.origin.r) - body.radius < max(self.width, self.height) / scale

    def set_scale(self, scale):
        """scale is only an indication. Compares it with maximum scale of visible bodies"""
        scales = [sprite.maxscale for sprite in self.spritegroup if self.is_visible(sprite.body, scale)]
        if scales:
            maxscale = min(scales)
            if 0 <= maxscale < scale:
                scale = maxscale
        if scale != self.prevscale:
            ScreenTransform.set_scale(self, scale)
            for i in self.spritegroup:
                if self.is_visible(i.body, scale):
                    i.rescale(scale)
            self.prevscale = scale
        return scale

    def add_body(self, body, imagename=None):
        r = re.compile(body.name + r"([0-9]*).(jpg|png)", re.I)
        if imagename is None:
            for i in os.listdir('images'):
                a = r.match(i)
                if a:
                    filename = i
                    idx = a.group(1)
                    if idx:
                        actradius = int(idx)
                    else:
                        actradius = -1
                    break
            else:
                filename = 'callisto.jpg'
                actradius = -1
        sprite = BodySprite(body, self, filename, actradius)
        self.spritegroup.add(sprite)
        sprite.rescale(self.scale)

    def add_spaceship(self, spaceship, imagelist=None):
        if imagelist is not None:
            imageidle, imageforward, imagebackward, imageexplosion = imagelist
        else:
            imageidle, imageforward, imagebackward, imageexplosion = 'rocket.png', 'rocket_lit.png', 'rocket_back.png', 'explosion.png'
        self.spaceshipsprite = SpaceShipSprite(spaceship, self, imageidle, imageforward, imagebackward, imageexplosion)
        self.spritegroup.add(self.spaceshipsprite)
        self.spaceshipsprite.rescale(self.scale)

    def add_bodylist(self, bodylist):
        for i in bodylist:
            if isinstance(i, cnbody.SpaceShip):
                self.add_spaceship(i)
            else:
                self.add_body(i)

    def clear(self):
        self.spritegroup.clear(self.screensurface, self.backgroundimage)

    def draw(self):
        self.spritegroup.draw(self.screensurface)

    def update(self):
        self.spritegroup.update()

    def onceinawhile(self):
        pass


class Text(pygame.sprite.DirtySprite):
    def __init__(self, text, font, color, width, height):
        # Call the parent class (Sprite) constructor  
        pygame.sprite.DirtySprite.__init__(self)

        self.font = font
        self.image = self.font.render(text, 1, color)
        self.Surf = pygame.Surface((width, height))
        W = self.image.get_width()
        H = self.image.get_height()
        self.Surf.blit(self.image, (height / 2 - H / 2, width / 2 - W / 2))
        self.visible = False
        self.rect = self.Surf.get_rect()

    def set_pos(self, x, y):
        self.rect.center = x, y


class BodyLabel(Text):
    def __init__(self, body, font):
        Text.__init__(self, body.name, font, (255, 255, 255), 100, 20)
        self.body = body


class Radar(Frame):

    def __init__(self, *args, **kwargs):
        self.font = loader.load_font(FONT_NAME, 10)
        self.textspritelist = []
        Frame.__init__(self, *args, **kwargs)

    def convert(self, pos):
        l, angle = cmath.polar(pos - self.origin.r)
        angle = angle - self.spaceshipsprite.body_angle / 180. * math.pi
        pos = cmath.rect(self.radar_formula(l), angle)
        return (self.sx_position + int(pos.real),
                self.height - self.sy_origin - int(pos.imag))

    def radar_formula(self, l):
        # return max(0,math.log(max(1,l))-14)
        return .4 * max(1, l) ** .2

    def add_spaceship(self, spaceship, imagelist=None):
        super(Radar, self).add_spaceship(spaceship, imagelist)
        self.spaceshipsprite.rotate = False

    def add_body(self, body, *args):
        Frame.add_body(self, body, *args)
        s = BodyLabel(body, self.font)
        s.visible = False
        self.spritegroup.add(s)
        self.textspritelist.append(s)

    def show_names(self):
        pos = {}
        for sprite in self.textspritelist:
            sprite.visible = self.safe_set_pos(sprite)
            x, y = sprite.rect.center
            if (x / 4, y / 4) in pos:
                sprite.visible = False
            else:
                pos[x / 4, y / 4] = True
            sprite.dirty = True

    def hide_names(self):
        for sprite in self.textspritelist:
            sprite.visible = False
            sprite.dirty = True


def log_quantize(x, q):
    lq = math.log(q)
    return math.exp(lq * int(math.log(x) / lq))


class LayeredFrame(Frame):
    def change_layer(self, sprite, layer):
        self.spritegroup.change_layer(sprite, layer)


class AutoZoomFrame(Frame):

    def __init__(self, name, screensurface, origin, autozoom, viewradius, backgroundimage, autozoom_timestep_mult,
                 backgroundcolor=None, gameloop=None):
        self.autozoom = False
        self.autozoom_timestep_mult = autozoom_timestep_mult
        Frame.__init__(self, name, screensurface, origin, viewradius, backgroundimage, backgroundcolor,
                       spritegroup=pygame.sprite.LayeredDirty())
        self.targetscale = self.scale
        self.currentscale = self.scale
        self.startmovescale = 1
        self.gameloop = gameloop
        self.set_autozoom(autozoom)

    def set_autozoom(self, autozoom):
        self.autozoom = autozoom
        if not autozoom and self.gameloop and self.gameloop.spaceship:
            self.origin = self.gameloop.spaceship

    def set_scale(self, scale):
        if self.autozoom:
            if scale != self.targetscale:
                self.targetscale = scale
                self.startmovescale = self.scale
        else:
            Frame.set_scale(self, scale)

    def update(self):
        if self.autozoom:
            if self.targetscale != self.scale:
                currentscale = self.scale + (self.targetscale - self.scale) / 2
                if abs(self.targetscale - currentscale) < abs(self.targetscale - self.startmovescale) * 0.2:
                    currentscale = self.targetscale
                Frame.set_scale(self, currentscale)
            # check for visibility
            if hasattr(self, 'spaceshipsprite') and not self.spaceshipsprite.visible:
                self.adapt_viewradius()
        Frame.update(self)

    def onceinawhile(self):
        """check sometimes if the zoom needs to be adapted"""
        if self.autozoom and self.gameloop:
            self.adapt_viewradius()

    def adapt_viewradius(self):
        """adapts screen resolution and simulation speed to location of spaceship compared
        to its closest body"""
        gameloop = self.gameloop
        if gameloop.spaceship:
            spaceship_r = gameloop.spaceship.r
            # finds where to place the origin of the view
            closest = i = gameloop.system.closest_body(gameloop.spaceship)
            # if the next parent is not too far, use it instead
            distance = closest_distance = abs(spaceship_r - closest.r) - closest.radius
            while True:
                if i.parent is None:
                    break
                new_distance = abs(spaceship_r - i.parent.r) - i.parent.radius
                if new_distance > closest_distance * 4:
                    break
                i = i.parent
                distance = new_distance
            # if the body at the origin is too large, use the projection of the spaceship on the surface
            if distance < i.radius:
                i = OnSurface(i, gameloop.spaceship)
                distance = abs(spaceship_r - i.r)
            # quantize the radius (otherwise would change permanently)
            targetradius = 1.2 * AUTOZOOM_QUANTIZE * max(log_quantize(distance, AUTOZOOM_QUANTIZE), closest.radius * .1)
            self.origin = i
            self.set_viewradius(targetradius)
            # magic formula to adapt the simulation speed according to the scale: accelerate when the
            # spaceship is far from any other body
            self.gameloop.set_time_step(
                self.autozoom_timestep_mult * max(1, 0.01 * log_quantize(self.targetscale ** -0.5, 2)))


class CollisionException(Exception):
    def __init__(self, bodylist):
        self.bodylist = bodylist

    def __str__(self):
        return 'Collision between\n' + '\n'.join(
            "%s and %s at %.2f km/h" % (i.name, j.name, (abs(i.v - j.v) * 3.6)) for i, j in self.bodylist)


CHEAT_SHEET = """
a	Autozoom
b/v	Inc/dec boot
[/]	Zoom in/out
esc	Quit
left/right	Turn 
up	Forward 
down	Reverse
p	Pause
l	Faster simul
k	Slower simul
0-9,bckspce	Edit DeltaV
Return	Apply DeltaV
"""
CHEAT_SHEET = [i.split('\t') for i in CHEAT_SHEET.strip().splitlines()]


class SpaceShipGame(gl.GameLoop):

    def __init__(self, caption, origin_time=None):
        self.paused = False
        self.deltav_string = ''
        self.deltav_command = 0
        gl.GameLoop.__init__(self, caption, frame_rate=20)
        self.simulation_time = 0
        self.set_time_step(TIME_STEP)
        self.user_command_acc = 0.
        self.angle_command = 0.
        self.body_angle = 0.
        self.user_command_angle = 0.
        self.acc_mult = ACC_MULT0
        self.origin_time = origin_time or time.time()

        self._bind_keys()

    def load_system_file(self, filename):
        builder = cnbody.NamedSystemBuilder()
        if filename.endswith('.py'):
            builder.execfile(filename)
        else:
            builder.readfile(filename)
        system = GameNBodySystem(builder.bodylist, COLL_TIME_MULTIPLIER)
        self.setsystem(system)

    def build_displays(self):
        """initializes the views"""
        backgroundimage = loader.load_image("background.jpg")
        radarimage = loader.load_image("radar.jpg")

        screen_rect = self.screensurface.get_rect()
        screen_height, screen_width = screen_rect.height, screen_rect.width
        radar_height = dashboard_width = 300
        cheatsheet_height = 200
        dashboard_height = screen_height - radar_height - cheatsheet_height
        if dashboard_height < 300:
            cheatsheet_height = 0
            dashboard_height = screen_height - radar_height

        mainframesurface = self.screensurface.subsurface(
            pygame.Rect(dashboard_width, 0, screen_width - dashboard_width, screen_height))
        origin = self.spaceship or Barycenter(self.system)
        max_radius = numpy.max(numpy.absolute(self.system.r - origin.r))
        self.mainframe = AutoZoomFrame("mainscreen",
                                       mainframesurface,
                                       origin,
                                       True,
                                       max_radius,
                                       backgroundimage, AUTOZOOM_TIMESTEP_MULT, gameloop=self)
        radarsurface = self.screensurface.subsurface(pygame.Rect(0, dashboard_height, dashboard_width, radar_height))
        self.radarframe = Radar("radar", radarsurface, DEFAULT_POS, 200e9, radarimage)

        display_list = [("time", 16), ("time acc", 10), ("auto zoom", 10)]
        if self.spaceship:
            display_list += [("closest", 16),
                             ("rel speed", 16), ("acc", 16),
                             ("radial v", 10), ("thrust", 13),
                             ("fuel", 16), ("gravity", 16),
                             ("attractor", 16),
                             ("escape", 13), ("eccentricity", 13),
                             ("delta v", 18)]
        self.display = DashboardDisplay(pygame.rect.Rect(0, 0, dashboard_width, dashboard_height), display_list)
        displayframesurface = self.screensurface.subsurface(pygame.Rect(0, 0, dashboard_width, dashboard_height))
        displayframe = SimpleFrame("display", displayframesurface, backgroundimage)
        displayframe.add(self.display)

        self.mainframe.add_bodylist(self.system.bodylist)
        self.radarframe.add_bodylist(self.system.bodylist)

        if self.spaceship:
            # self.mainframe.origin=self.system.closest_body(self.spaceship)
            self.radarframe.origin = self.spaceship
            self.add_frames(self.mainframe, self.radarframe, displayframe)
        else:
            #    for i in self.system.bodylist:
            #        if i.parent is None:
            #            self.mainframe.origin = i
            self.add_frames(displayframe, self.mainframe)

        if cheatsheet_height > 0:
            cheatsheet = DashboardDisplay(pygame.rect.Rect(10, 0, dashboard_width - 10, cheatsheet_height),
                                          [(i[0], 10) for i in CHEAT_SHEET],
                                          refreshfreq=200,
                                          margin=1)
            for k, v in CHEAT_SHEET:
                cheatsheet.values[k] = v
            cheatsheetframesurface = self.screensurface.subsurface(
                pygame.Rect(0, dashboard_height + radar_height, dashboard_width, cheatsheet_height))
            cheatsheetframe = SimpleFrame("cheatsheet", cheatsheetframesurface, backgroundimage)
            cheatsheetframe.add(cheatsheet)
            self.add_frames(cheatsheetframe)

        if self.mainframe.autozoom:
            self.mainframe.adapt_viewradius()

    def _bind_keys(self):
        self.bind(QUIT, self.quit)

        self.bindkey(KEYDOWN, K_ESCAPE, self.quit)
        self.bindkey(KEYDOWN, K_RIGHT, self.on_rotateright)
        self.bindkey(KEYDOWN, K_LEFT, self.on_rotateleft)
        self.bindkey(KEYDOWN, K_DOWN, self.on_brake)
        self.bindkey(KEYDOWN, K_UP, self.on_accelerate)
        self.bindkey(KEYDOWN, K_F11, self.toggle_fullscreen)
        self.bindkey(KEYUP, K_UP, self.on_stop)
        self.bindkey(KEYUP, K_DOWN, self.on_stop)
        self.bindkey(KEYUP, K_RIGHT, self.on_stoprotation)
        self.bindkey(KEYUP, K_LEFT, self.on_stoprotation)

        self.bindkey(KEYDOWN, K_p, self.on_pause)
        self.bindkey(KEYDOWN, K_a, self.toggle_autozoom)
        self.bindkey(KEYDOWN, K_l, self.faster)
        self.bindkey(KEYDOWN, K_k, self.slower)
        self.bindkey(KEYDOWN, K_b, self.on_increasethrust)
        self.bindkey(KEYDOWN, K_v, self.on_reducethrust)
        self.bindkey(KEYPRESSED_EVENT, K_LEFTBRACKET, self.zoom_in)
        self.bindkey(KEYPRESSED_EVENT, K_RIGHTBRACKET, self.zoom_out)
        self.bindfigures(self.set_deltav)
        self.bindkey(KEYDOWN, K_PERIOD, self.set_deltav)
        self.bindkey(KEYDOWN, K_KP_PERIOD, self.set_deltav)
        self.bindkey(KEYDOWN, K_BACKSPACE, self.del_deltav)
        self.bindkey(KEYPRESSED_EVENT, K_RETURN, self.apply_deltav)
        self.bindkey(KEYPRESSED_EVENT, K_KP_ENTER, self.apply_deltav)
        self.bindkey(KEYDOWN, K_n, self.on_show_names)
        self.bindkey(KEYUP, K_n, self.on_hide_names)

    def set_time_step(self, time_step):
        if not self.paused:
            self.time_step = time_step

    def set_deltav(self, event):
        dv = self.deltav_string + chr(event.key)
        if dv != '.':
            try:
                float(dv)
            except:
                return
        self.deltav_string = dv

    def del_deltav(self, event):
        self.deltav_string = self.deltav_string[:-1]

    def apply_deltav(self, event):
        if self.deltav_string:
            self.deltav_command = float(self.deltav_string)
            self.deltav_string = ''

    def on_show_names(self, event):
        self.radarframe.show_names()

    def on_hide_names(self, event):
        self.radarframe.hide_names()

    def on_accelerate(self, event):
        self.user_command_acc = self.acc_mult

    def on_brake(self, event):
        self.user_command_acc = -self.acc_mult

    def on_increasethrust(self, event):
        if self.acc_mult < 40:
            self.acc_mult *= 1.5

    def on_reducethrust(self, event):
        if self.acc_mult > 0.1:
            self.acc_mult /= 1.5

    def on_stop(self, event):
        self.user_command_acc = 0.

    def on_pause(self, event):
        if self.paused:
            self.time_step = self.original_time_step
        else:
            self.original_time_step = self.time_step
            self.time_step = 0
        self.paused = not self.paused

    def on_rotateleft(self, event):
        self.user_command_angle = ROTATION_DEG_BY_SEC

    def on_rotateright(self, event):
        self.user_command_angle = -ROTATION_DEG_BY_SEC

    def on_stoprotation(self, event):
        self.user_command_angle = 0.

    def toggle_autozoom(self, event):
        self.mainframe.set_autozoom(not self.mainframe.autozoom)

    def zoom_in(self, event):
        self.mainframe.set_scale(self.mainframe.scale * 1.1)

    def zoom_out(self, event):
        self.mainframe.set_scale(self.mainframe.scale / 1.1)

    def faster(self, event):
        if self.mainframe.autozoom:
            self.mainframe.autozoom_timestep_mult *= 1.5
        else:
            self.set_time_step(self.time_step * 2)

    def slower(self, event):
        if self.mainframe.autozoom:
            self.mainframe.autozoom_timestep_mult /= 1.5
        else:
            self.set_time_step(self.time_step / 2)

    def setsystem(self, system):
        self.system = system
        self.dt = self.recommanded_dt()

        # find the spaceship
        self.spaceship = None
        for i in system.bodylist:
            if i.name == 'Spaceship':
                self.spaceship = i
                self.spaceship.fuel = FUEL

    def recommanded_dt(self):
        """recommanded dt for the simulation - this is not the time step"""
        try:
            return 0.01 * self.system.collision_time()
        except ValueError:
            return 1

    def pre_stepforward(self):
        if self.paused:
            return
        self.simul_time_step = self.clock.get_time() * 0.001  # time since last frame, in s
        self.dt = self.recommanded_dt()
        self.body_angle += min(self.user_command_angle * self.simul_time_step, 15)
        if self.spaceship:
            if self.deltav_command > 0:
                correct_acc = min(40, self.deltav_command / self.time_step)
                self.deltav_command -= self.time_step * correct_acc
            else:
                correct_acc = self.user_command_acc  # *min(1,400/self.get_time_acc())
            self.spaceship.acc_command_polar = correct_acc, self.body_angle
            self.spaceship.fuel -= abs(correct_acc * self.time_step)

            for i in self.framelist:
                if hasattr(i, 'spaceshipsprite'):
                    i.spaceshipsprite.user_command_acc = correct_acc
                    i.spaceshipsprite.body_angle = self.body_angle

    def step(self):
        try:
            self.pre_stepforward()

            if self.time_step > 0:
                steps, self.simulation_time = self.system.forward_fixed(self.simulation_time,
                                                                        self.simulation_time + self.time_step,
                                                                        self.dt)
                self.post_stepforward(self.time_step)
        except CrashException:
            if self.spaceship:
                self.mainframe.origin = self.spaceship.closest
                self.mainframe.set_viewradius(abs(self.spaceship.r - self.spaceship.closest.r) * 2)
                self.mainframe.spaceshipsprite.crash()
                self.spaceship = None

    def get_time_acc(self):
        """factor of time acceleration"""
        return max(1e-6, self.time_step * 1000. / self.clock.get_time())

    def post_stepforward(self, dt):

        self.display.values['time'] = time.strftime("%x %X", time.localtime(self.origin_time + self.simulation_time))
        self.display.values['time acc'] = "%.1f" % self.get_time_acc()
        self.display.values['auto zoom'] = "on" if self.mainframe.autozoom else "off"

        if self.spaceship:
            closest = self.spaceship.closest
            attractor = self.spaceship.attractor
            rel_velocity = self.spaceship.v - closest.v
            rel_pos = self.spaceship.r - closest.r
            self.display.values['closest'] = "%s   %s" % (closest.name, format_distance(abs(rel_pos) - closest.radius))
            self.display.values['rel speed'] = format_speed(abs(rel_velocity))
            self.display.values['acc'] = '%.2f' % (abs(self.spaceship.a))
            self.display.values['radial v'] = format_speed(self.spaceship.rel_coord.v.real)
            self.display.values['thrust'] = '%.2f' % (self.acc_mult)
            self.display.values['fuel'] = '%.2f' % (self.spaceship.fuel)
            self.display.values['gravity'] = '%.2f' % (abs(self.spaceship.a - self.spaceship.acc_command))
            attractor_string = attractor.name
            if attractor != closest:
                attractor_string += '  ' + format_distance(abs(self.spaceship.r - attractor.r) - attractor.radius)
            self.display.values['attractor'] = attractor_string
            self.display.values['escape'] = '%.1f' % self.spaceship.escape_speed(attractor)
            ecc = self.spaceship.eccentricity(attractor)
            if ecc > 10:
                self.display.values['eccentricity'] = '--'
            else:
                self.display.values['eccentricity'] = "%.2f" % ecc
            self.display.values['delta v'] = self.deltav_string or str(self.deltav_command)


def main(filename):
    pygame.init()
    gameloop = SpaceShipGame("Orbital Flight")

    # read the system from file
    gameloop.load_system_file(filename)
    gameloop.build_displays()
    pygame.display.update()
    # and loop for ever or until something happens
    gameloop.mainloop()


# os.environ['SDL_VIDEODRIVER'] = 'windib'

# this calls the 'main' function when this script is executed
if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main('data/planet.csv')
