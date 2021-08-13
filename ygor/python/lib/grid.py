import pygame
import random

import lib.colors as colors

from lib.vector import Vector
from lib.utils import bilinear_interpolation, bicubic_interpolation

class Cell:

    def __init__(self, p, coord, block_size):
        self.p = p
        self.coord = coord
        self.size = block_size
        l = list(range(-10, 11))
        l.remove(0)
        u, v = random.choice(l), random.choice(l)
        self.vectors = [
            Vector(self.center, self.size//4, u, v),
            Vector((self.x, self.y+self.size//2), self.size//4, u, 0, 'Staggered'),
            Vector((self.x+self.size//2, self.y+self.size), self.size//4, 0, v, 'Staggered')
        ]
    
    @property
    def x(self):
        return self.coord[0]

    @property
    def y(self):
        return self.coord[1]
    
    @property
    def center(self):
        return (self.x+self.size//2, self.y+self.size//2)
    
    def draw_vectors(self, screen, mode):
        for v in self.vectors:
            if v.mode == mode:
                v.draw(screen)

class Grid:

    def __init__(self, width, height, block_size):
        self.width = width
        self.height = height
        self.blocks = []
    
        self.generate_grid(block_size)
    
    def generate_grid(self, block_size):
        self.blocks = [[Cell(random.randint(1, 10), (x, y), block_size) 
            for y in range(0, self.height, block_size)] 
            for x in range(0, self.width, block_size)]
    
    def draw_grid(self, screen):
        for row in self.blocks:
            for cell in row:
                rect = pygame.Rect(cell.x, cell.y, cell.size, cell.size)
                color = (0,0, cell.p * 25)
                screen.fill(color, rect)
                pygame.draw.rect(screen, colors.BLACK, rect, 1)

    def draw_vectors(self, screen, mode):
        for row in self.blocks:
            for cell in row:
                cell.draw_vectors(screen, mode)
    
    def add_vector(self, pos, interp_mode):
        b_size = self.blocks[0][0].size

        i, j = pos[0]//b_size, pos[1]//b_size
        block = self.blocks[i][j]

        if interp_mode == 'Bilinear':
            u, v = bilinear_interpolation(i, j, pos, self.blocks)
        elif interp_mode == 'Bicubic':
            u, v = bicubic_interpolation(i, j, pos, self.blocks)
        
        block.vectors.append(Vector(pos, block.size//4, u, v))
        block.vectors.append(Vector(pos, block.size//4, u, 0, 'Staggered'))
        block.vectors.append(Vector(pos, block.size//4, 0, v, 'Staggered'))
