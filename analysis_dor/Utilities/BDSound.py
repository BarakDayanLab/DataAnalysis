import os
from enum import Enum, auto
from playsound import playsound


class SOUNDS(Enum):
    DEFAULT = 0  # auto()
    INDICATOR_1 = 1  # auto()
    INDICATOR_2 = 2  # auto()
    WARN = 3  # auto()
    FATAL = 4  # auto()


class BDSound:

    def __init__(self):
        """
        The class will initiate the media folder
        """
        # TODO: make this dependant on the user - are we on what machine?
        self.media_folder = r'C:\Windows\Media'
        pass

    def play(self, sound_enum):

        name = 'ding'  # default sound
        if sound_enum == SOUNDS.INDICATOR_1:
            name = 'ding'
        elif sound_enum == SOUNDS.INDICATOR_2:
            name = 'chord'
        elif sound_enum == SOUNDS.WARN:
            name = 'chimes'
        elif sound_enum == SOUNDS.FATAL:
            name = 'Alarm05'
        else:
            print(f'Unknown sound type: {sound_enum}')

        sound_file = os.path.join(self.media_folder, name + '.wav')
        playsound(sound_file)
        pass
