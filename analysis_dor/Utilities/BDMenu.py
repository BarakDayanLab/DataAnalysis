import os
import json
import sys
import numpy as np


class BDMenu:

    def __init__(self, caller, menu_file, menu_json=None):

        self.caller = caller

        if menu_file:
            self._load_menu(menu_file)
        elif menu_json:
            self.menu_data = menu_json

        pass

    def _load_menu(self, menu_file):

        try:
            # Opening JSON file
            f = open(menu_file)

            # returns JSON object as
            # a dictionary
            self.menu_data = json.load(f)

            # Closing file
            f.close()
        except Exception as err:
            print(err)

    def func_not_found(self):  # just in case the menu code cannot invoke a function in the caller object
        print ('No Function Found!')

    def display(self):

        while True:
            self._display()
            if "infinite" not in self.menu_data or self.menu_data['infinite']==False:
                break
        print('Bye!')

    def _display(self):

        # TODO: check validity of json menu items

        call_map = {}

        # Iterating through the json menu items
        for menu_item in self.menu_data['menu_items']:
            if 'spacer' in menu_item:
                print()
            else:
                print(f'{menu_item["order"]}) {menu_item["display"]}')
                call_map[menu_item["order"]] = menu_item

        # Wait for input
        selection = input(">> ")
        selection = int(selection)

        # Sanity check
        if "action" not in call_map[selection]:
            self._error('Missing "action" section in menu definitions. Check your json file.')

        # Get the function name and invoke it
        func_name = call_map[selection]["action"]

        # If it is an exit function, then exit
        if func_name.lower() == 'exit':
            os._exit(0)
            return

        args = {}
        if "args" in call_map[selection]:
            for arg in call_map[selection]["args"]:
                prompt = arg['name']
                # Is there a default value suggested?
                if arg['default']:
                    default_values_str = self._values_to_str(arg['default'])
                    prompt = f'{prompt} (Press <Enter> for Default: {default_values_str}) >>'
                value_str = input(prompt)
                # Should we be using default values? (if user just pressed <enter> w/o entering a value
                if len(value_str) == 0:
                    value_str = default_values_str
                value = self._convert_input(value_str, arg['type'])

                # Note - we always use lowercase for the arguments - so make sure the update function has lowercase args
                args[arg['name'].lower()] = value

        # Invoke the method/action
        func = getattr(self.caller, func_name, self.func_not_found)
        result = None
        try:
            if len(args) == 0:
                result = func()
            else:
                result = func(**args)
        except Exception as err:
            print(f"Unable to invoke function {func} or failed during invocation - {err}")

        return result

    def _values_to_str(self, values):
        return str(values)

    def _convert_input(self, value_str, value_type):
        # If value_str is encapsulated by brackets, we remove them
        if value_str.startswith('[') and value_str.endswith(']'):
            value_str = value_str[1:-1]
        if 'array_of' in value_type:
            # Remove redundant spaces and split by commas
            value_str = value_str.replace(" ", "").split(',')
            # Get the type requested
            value_type = value_type.replace('array_of_', '')
            # Create array
            values = []
            for val in value_str:
                values.append(self._convert_input(val, value_type))
            return values
        elif value_type == 'int' or value_type == 'integer':
            return int(float(value_str))  # Python fails to convert float strings with this: int("3.1415"), so we need to convert to floar first and then int
        elif value_type == 'float':
            return float(value_str)
        elif value_type == 'string' or value_type == 'str':
            return str(value_str)
        else:
            self._error(f'Unknown value_type ({value_type}). Check your json.')
        pass
    def _error(self, message):
        print(message)
        sys.exit("User exited")
        pass

    def unit_test_callback(self, float_values):
        pass

    def unit_tests(self):
        menu_json = {
            "menu_items": [
                {
                    "order": 1,
                    "display": "Measure Temperature",
                    "action": "unit_test_callback",
                    "args": [
                        {
                            "name": "float_values",
                            "type": "array_of_float",
                            "default": "[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.5]"
                        },
                    ]
                },
                {
                    "order": 9,
                    "display": "Exit",
                    "action": "exit"
                }
            ]
        }
        self.float_values = [2.0, 4.0]
        menu = BDMenu(caller=self, menu_file=None, menu_json=menu_json)
        menu.display()
        pass

if __name__ == "__main__":

    BDMenu(caller=None, menu_file=None, menu_json=None).unit_tests()
