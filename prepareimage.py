import pygame

class Image(object):
    def __init__(self,name,size=None):
        self.image = pygame.image.load(name)
        if not size:
            size = min(self.image.get_rect().size)
        radius = size/2
        center = self.image.get_rect().center



