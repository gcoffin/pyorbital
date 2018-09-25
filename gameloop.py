import pygame
from pygame.locals import *
ONCE_IN_A_WHILE=10
KEYPRESSED_EVENT=42

SCREEN_WIDTH, SCREEN_HEIGHT = (1200, 800)

'''a simple game framework that provides some support
to manage events and the loop'''

class eventhandler(object):
    def __init__(self):
        self.map={}
        self.pressedkeys=set()

    def bind(self,eventtype,action):
        self.map.setdefault(eventtype,[]).append(action)

    def bindkey(self,eventtype,key,action):
        if eventtype not in self.map:
            self.map[eventtype]=keyeventhandler()
        self.map[eventtype].addkey(key,action)

    def bindchars(self,chars,action):
        '''binds an action to specified chars'''
        for c in chars:
            self.bindkey(KEYDOWN,globals()['K_'+c],action)

    def bindfigures(self,action):
        '''binds an action to keys 0-9'''
        self.bindchars([str(i) for i in range(10)],action)
        self.bindchars(['KP%d'%i for i in range(10)],action)

    def __call__(self,event):
        for action in self.map.getdefault(event,[]):
            action(event)

    def pumpevents(self):
        '''calls registerd event handlers according to events in the
        queue'''
        for event in pygame.event.get():
            if event.type==KEYDOWN:
                self.pressedkeys.add(event.key)
            if event.type==KEYUP:
                self.pressedkeys.discard(event.key)
            for action in self.map.get(event.type,[]):
                action(event)
        # special management for pressed keys
        for key in self.pressedkeys:
            for action in self.map.get(KEYPRESSED_EVENT,[]):
                event=pygame.event.Event(KEYPRESSED_EVENT,key=key)
                action(event)


class keyeventhandler(object):
    def __init__(self):
        self.map={}

    def addkey(self,key,action):
        self.map.setdefault(key,[]).append(action)

    def __call__(self,event):
        if event.key in self.map:
            for action in self.map[event.key]:
                action(event)

    def __iter__(self):
        return iter([self])

class GameLoop(eventhandler):
    def __init__(self,caption,frame_rate):
        eventhandler.__init__(self)
        self.frame_rate = frame_rate
        self.fullscreen = True
        self.caption = caption
        self.set_display_mode()
        self.commandbind=keyeventhandler()
        self.framelist=[]
        self.paused = False

    def add_frames(self,*frames):
        self.framelist.extend(frames)

    def pre_stepforward(self):
        pass

    def post_stepforward(self,dt):
        pass

    def onceinawhile(self):
        pass

    def quit(self,event):
        self.running=False

    def mainloop(self):
        t=0
        # Initialize some control variables
        self.running = True
        self.clock = pygame.time.Clock()
        countkey=0

        pygame.display.flip()

        for f in self.framelist:
            f.onceinawhile()
        onceinawhilecount=0

        while self.running is True:
            self.clock.tick( self.frame_rate )
            for frame in self.framelist:
                frame.clear()
            self.pumpevents()

            onceinawhilecount+=1

            if onceinawhilecount>=ONCE_IN_A_WHILE:
                self.onceinawhile()

                for f in self.framelist:
                    f.onceinawhile()

                onceinawhilecount=0

            self.step()

            for f in self.framelist:
                f.update()

            for f in self.framelist:
                f.draw()

            pygame.display.update()

    def set_display_mode(self):
        if self.fullscreen:
            self.screensurface = pygame.display.set_mode((0, 0), FULLSCREEN)
        else:
            self.screensurface = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT), 0, 32)
        pygame.display.set_caption(self.caption)

    def toggle_fullscreen(self,event=None):
        # todo: need to recalculate layout on toggle
        self.fullscreen = not self.fullscreen
        self.set_display_mode()
