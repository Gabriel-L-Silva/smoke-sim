import pygame
import pygame_gui

from pygame_gui.elements import UIWindow
from pygame_gui.elements import UILabel
from pygame_gui.elements import UIDropDownMenu

class Menu(UIWindow):

    def __init__(self, ui_manager):
        super().__init__(pygame.Rect((50, 50), (350, 200)), ui_manager,
                         window_display_title='Config',
                         object_id='#config')
    
        self.vector_mode_label = UILabel(pygame.Rect((5,5), (100, 25)),
                                        'Vetor Mode: ',
                                        self.ui_manager,
                                        container=self)
        
        self.current_vector_mode = 'Central'
        self.vector_mode = UIDropDownMenu(['Central', 'Staggered'],
                                            self.current_vector_mode,
                                            pygame.Rect((140, 5),
                                                        (100, 25)),
                                            self.ui_manager,
                                            container=self)
        
        self.interpolation_mode_label = UILabel(pygame.Rect((-8, 50), (150, 25)),
                                                    'Interpolation: ',
                                                    self.ui_manager,
                                                    container=self)

        self.current_interp_mode = 'Bilinear'
        self.interpolation_mode = UIDropDownMenu(['Bilinear', 'Bicubic'],
                                                self.current_interp_mode,
                                                pygame.Rect((140, 50),
                                                            (100, 25)),
                                                self.ui_manager,
                                                container=self)
        
        self.set_blocking(True)
        self.close_window_button.kill()

    def process_event(self, event: pygame.event.Event) -> bool:
        super().process_event(event)

        if event.type == pygame.USEREVENT:
            if event.user_type == pygame_gui.UI_DROP_DOWN_MENU_CHANGED:
                if event.ui_element == self.vector_mode:
                    self.current_vector_mode = self.vector_mode.selected_option
                elif event.ui_element == self.interpolation_mode:
                    self.current_interp_mode = self.interpolation_mode.selected_option