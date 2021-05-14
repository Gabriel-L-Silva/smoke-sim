import pygame
import pygame_gui
import sys

from lib.grid import Grid
from lib.menu import Menu

WINDOW_HEIGHT = 800
WINDOW_WIDTH = 1200
BLOCKSIZE = 100

def main():
    pygame.init()
    screen = pygame.display.set_mode((WINDOW_WIDTH, WINDOW_HEIGHT))
    clock = pygame.time.Clock()
    manager = pygame_gui.UIManager((WINDOW_WIDTH, WINDOW_HEIGHT))

    grid = Grid(WINDOW_WIDTH, WINDOW_HEIGHT, BLOCKSIZE)
    menu = Menu(manager)
    menu.hide()

    while True:
        time_delta = clock.tick(30)/1000

        grid.draw_grid(screen)
        grid.draw_vectors(screen, menu.current_vector_mode)

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()

            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_r:
                    grid.generate_grid(BLOCKSIZE)

            if event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == 1:  # left click
                    if not menu.visible:
                        pos = pygame.mouse.get_pos()
                        grid.add_vector(pos, menu.current_interp_mode)
                
                if event.button == 3:  # right click
                    if menu.visible:
                        menu.hide()
                    else:
                        menu.show()

            manager.process_events(event)
        
        manager.update(time_delta)
        manager.draw_ui(screen)

        pygame.display.update()


if __name__ == '__main__':
    main()