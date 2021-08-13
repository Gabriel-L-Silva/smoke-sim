import pygame
import math

from pygame import color

import lib.colors as colors

from lib.utils import draw_arrow

class Vector:

    def __init__(self, coord, length, u, v, mode='Central'):
        self.u = u
        self.v = v
        self.length = length
        self.coord = coord
        self.mode = mode

    @property
    def x(self):
        return self.coord[0]

    @property
    def y(self):
        return self.coord[1]    

    def draw(self, screen):
        x2, y2 = (int(self.x + self.length * (-1 if self.u < 0 else 1) * (1 if self.u else 0)), 
                  int(self.y + self.length * (-1 if self.v < 0 else 1) * (1 if self.v else 0)))
        
        thickness = (abs(self.u) + abs(self.v) + 1)//2
        
        draw_arrow(screen, self.coord, (x2, y2), thickness, colors.WHITE)
