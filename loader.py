import os, pygame
from pygame.locals import *
import pygame.font

DEFAULT_FILE = 'callisto.jpg'

DIRECTORY = os.path.dirname(__file__)

all_images={}

def load_image(name, colorkey=None, path=None):
    global all_images
    if name in all_images:
        return all_images[name]
    fullname = os.path.join(path or DIRECTORY,'images', name)
    try:
        image = pygame.image.load(fullname)
    except pygame.error as message:
        if name!=DEFAULT_FILE:
            return load_image(DEFAULT_FILE)
        print('Cannot load image:', name)
        raise SystemExit(message)

    image = image.convert_alpha()

    if colorkey is not None:
        if colorkey is -1:
            colorkey = image.get_at((0,0))
        image.set_colorkey(colorkey, RLEACCEL) # set transparent color

    all_images[name]=image

    return image

def load_font(name, size,path=None):
    if name is None:
        fullname = None
        size =int(size*1.5) # default font is too small
    else:
        fullname = os.path.join(path or DIRECTORY,'fonts', name)
    try:
        font = pygame.font.Font(fullname, size)
    except pygame.error as message:
        print( 'Cannot load font:', os.path.abspath(fullname))
        raise SystemExit( message)

    return font
    
