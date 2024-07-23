import os
import time

# TODO:
#  1) implement logging into file
#  3) Debug - filtering

class BDLogger:
    HEADER = '\033[95m'
    OK_BLUE_COLOR = '\033[94m'
    OK_CYAN_COLOR = '\033[96m'
    INFO_COLOR = '\033[92m'  # Green
    WARNING_COLOR = '\033[93m'  # Read
    FAIL_COLOR = '\033[91m'
    DEBUG_COLOR = '\033[96m'  # Cyan
    END_COLOR = '\033[0m'
    BOLD_FORMAT = '\033[1m'
    UNDERLINE_FORMAT = '\033[4m'

    def __init__(self):
        self.log_file = None  # Not saving the log
        pass

    def __del__(self):
        if self.log_file:
            self.log_file.close()

    def turn_save_on(self, log_path):
        log_file = os.path.join(log_path, 'log.txt')
        self.log_file = open(log_file, "a+")  # append mode
        pass

    def blue(self, str=''):
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        time_stamped_msg = f'{time_str} - {str}'
        formatted_msg = self._resolve_format(f'{BDLogger.OK_BLUE_COLOR}{time_stamped_msg}{BDLogger.END_COLOR}')
        self._print(time_stamped_msg, formatted_msg)

    def error(self, str=''):
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        time_stamped_msg = f'{time_str} - Error: {str}'
        formatted_msg = self._resolve_format(f'{BDLogger.FAIL_COLOR}{time_stamped_msg}{BDLogger.END_COLOR}')
        self._print(time_stamped_msg, formatted_msg)

    def warn(self, str=''):
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        time_stamped_msg = f'{time_str} - Warning: {str}'
        formatted_msg = self._resolve_format(f'{BDLogger.WARNING_COLOR}{time_stamped_msg}{BDLogger.END_COLOR}')
        self._print(time_stamped_msg, formatted_msg)

    def info(self, str=''):
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        time_stamped_msg = f'{time_str} - Info: {str}'
        formatted_msg = self._resolve_format(f'{BDLogger.INFO_COLOR}{time_stamped_msg}{BDLogger.END_COLOR}')
        self._print(time_stamped_msg, formatted_msg)

    def debug(self, str=''):
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        time_stamped_msg = f'{time_str} - Debug: {str}'
        formatted_msg = self._resolve_format(f'{BDLogger.DEBUG_COLOR}{time_stamped_msg}{BDLogger.END_COLOR}')
        self._print(time_stamped_msg, formatted_msg)

    def _print(self, time_stamped_msg, formatted_msg):
        if self.log_file:
            self.log_file.write(time_stamped_msg + '\n')
            self.log_file.flush()
        print(formatted_msg)
        pass

    def _resolve_format(self, str):
        """ TODO: fix bug where we do not return the color format"""
        str = str.replace('**', BDLogger.BOLD_FORMAT, 1)
        str = str.replace('**', BDLogger.END_COLOR, 2)
        str = str.replace('__', BDLogger.UNDERLINE_FORMAT, 1)
        str = str.replace('__', BDLogger.END_COLOR, 2)
        return str