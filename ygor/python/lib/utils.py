import pygame
import math
import random


def bilinear_aux(x, y, points):
    # http://en.wikipedia.org/wiki/Bilinear_interpolation

    points = sorted(points)  # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError('points do not form a rectangle')
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        raise ValueError('(x, y) not within the rectangle')

    r = (q11 * (x2 - x) * (y2 - y) +
         q21 * (x - x1) * (y2 - y) +
         q12 * (x2 - x) * (y - y1) +
         q22 * (x - x1) * (y - y1)
        ) / ((x2 - x1) * (y2 - y1))
    
    return round(r)


def bilinear_interpolation(i, j, pos, blocks):
    points_u = ()
    points_v = ()
    size = blocks[0][0].size

    try:
        points_u = (
            (blocks[i][j].x, blocks[i][j].y, blocks[i][j].vectors[0].u),
            (blocks[i+1][j].x+size, blocks[i+1][j].y, blocks[i+1][j].vectors[0].u),
            (blocks[i][j+1].x, blocks[i][j+1].y+size, blocks[i][j+1].vectors[0].u),
            (blocks[i+1][j+1].x+size, blocks[i+1][j+1].y+size, blocks[i+1][j+1].vectors[0].u),
        )
        points_v = (
            (blocks[i][j].x, blocks[i][j].y, blocks[i][j].vectors[0].v),
            (blocks[i+1][j].x+size, blocks[i+1][j].y, blocks[i+1][j].vectors[0].v),
            (blocks[i][j+1].x, blocks[i][j+1].y+size, blocks[i][j+1].vectors[0].v),
            (blocks[i+1][j+1].x+size, blocks[i+1][j+1].y+size, blocks[i+1][j+1].vectors[0].v),
        )
    except IndexError:
        points_u = (
            (blocks[i-1][j-1].x, blocks[i-1][j-1].y, blocks[i-1][j-1].vectors[0].u),
            (blocks[i][j-1].x+size, blocks[i][j-1].y, blocks[i][j-1].vectors[0].u),
            (blocks[i-1][j].x, blocks[i-1][j].y+size, blocks[i-1][j].vectors[0].u),
            (blocks[i][j].x+size, blocks[i][j].y+size, blocks[i][j].vectors[0].u),
        )
        points_v = (
            (blocks[i-1][j-1].x, blocks[i-1][j-1].y, blocks[i-1][j-1].vectors[0].v),
            (blocks[i][j-1].x+size, blocks[i][j-1].y, blocks[i][j-1].vectors[0].v),
            (blocks[i-1][j].x, blocks[i-1][j].y+size, blocks[i-1][j].vectors[0].v),
            (blocks[i][j].x+size, blocks[i][j].y+size, blocks[i][j].vectors[0].v),
        )
    finally:
        return bilinear_aux(pos[0], pos[1], points_u), bilinear_aux(pos[0], pos[1], points_v)


def bicubic_interpolation(i, j, pos, blocks):
    l = list(range(-10,11))
    l.remove(0)
    return (random.choice(l), random.choice(l))


def draw_arrow(screen, start, end, thickness, color):
    rad = math.pi/180
    pygame.draw.line(screen, color, start, end, thickness)

    rotation = (math.atan2(start[1] - end[1], end[0] - start[0])) + math.pi/2
    arrow_thick = round(thickness*1.2)
    pygame.draw.polygon(screen, color, ((end[0] + arrow_thick * math.sin(rotation),
                                         end[1] + arrow_thick * math.cos(rotation)),
                                        (end[0] + arrow_thick * math.sin(rotation - 120*rad),
                                         end[1] + arrow_thick * math.cos(rotation - 120*rad)),
                                        (end[0] + arrow_thick * math.sin(rotation + 120*rad),
                                         end[1] + arrow_thick * math.cos(rotation + 120*rad))))
    