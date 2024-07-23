import sys
from Utilities.Utils import Utils


class BDKeyboard:

    def __init__(self, keyboard_mapping):

        self.listener = None

        # Set keyboard definitions
        self.keyboard_mapping = keyboard_mapping

        # Convert it to a map of entries by keystrokes
        self.callbacks = {}
        for attribute, entry in self.keyboard_mapping.items():
            for key in entry['keys']:
                self.callbacks[key] = entry
        pass

    def invoke_callback_func(self, event):
        func_name = ''
        # Find the corresponding key-press
        if event.key in self.callbacks:
            try:
                func_name = 'keyboard_handler__' + self.callbacks[event.key]['callback']
                func = getattr(self.listener, func_name)
                func(event.key)
            except AttributeError:
                print(f'Could not find function {func_name} in listener. Check your spelling.')
            except Exception as err:
                print(f'Invocation error in function {func_name}, {err}')
        pass

    def set_listener(self, figure, listener):
        self.listener = listener
        figure.canvas.mpl_connect('key_press_event', self.on_press)

    def on_press(self, event):
        self.invoke_callback_func(event)
        #print(f'pressed: [{event.key}]')
        sys.stdout.flush()