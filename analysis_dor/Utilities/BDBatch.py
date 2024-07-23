import os
import json
from Utilities.Utils import Utils


class BDBatch:

    def __init__(self, json_map_path= None, batches_map=None, batch_size=None):

        if batches_map is not None:
            self.batches_map = batches_map
        else:
            try:
                the_path = '.' if json_map_path is None else json_map_path
                the_file = os.path.join(the_path, 'batches_map.json')
                f = open(the_file)
                self.batches_map = json.load(f)
                f.close()
            except Exception as err:
                print(f'Unable to open batches_map.json. Reason: {err}')

        self.batch_size = batch_size

        # Initialize the empty batches
        self.empty_all()
        pass

    def __getitem__(self, key):
        if key not in self.batches_map:
            raise Exception(f'KeyError: no key {key} in batcher')
        return self.batches_map[key]['batched_data']

    def __setitem__(self, key, new_value):
        if key not in self.batches:
            self.batches_map[key]['batched_data'] = new_value

    def __delitem__(self, key):
        del self.batches_map[key]

    def empty_batch(self, key):
        # TODO: handle case where this is an array
        self.batches_map[key]['batches_data'] = []

    def set_batch_size(self, batch_size):
        self.batch_size = batch_size

    def empty_all(self):
        for item in self.batches_map.values():
            item['batched_data'] = []
        pass

    def add(self, key, value):
        if key not in self.batches_map:
            print(f'Cannot find batch {key} in batches. Not adding.')
        self.batches_map[key]['batched_data'].append(value)

    def batch_all(self, self_object):
        """
        This function scans the SELF object it got as an argument.
        It finds all members that are relevant to batch (as defined)
        Each member is batched into the relevant batch.
        """
        values_names = [x['value_name'] for x in self.batches_map.values() if 'value_name' in x]
        values_dict = Utils.get_self_members_and_values(self_object, specific=values_names)

        for batch_name, entry in self.batches_map.items():

            # If there's no value_name mentioned, we ignore this batch, it is probably done by a specific call
            if "value_name" not in entry:
                continue

            value_name = entry["value_name"]

            # Ensure the value is indeed in self object
            if value_name not in values_dict:
                #print(f'WARN: Check your batches map - no value {value_name} in self object. Not batching it.')
                continue

            # Get the relevant value we want to add to batch
            val = values_dict[value_name]

            # Check if it's empty and we're prevented from adding empty values
            if "allow_empty" in entry and not entry["allow_empty"] and self._is_value_empty(val) == 0:
                continue

            # Add the value (either FIFO style or at the end)
            if self.batch_size is None:
                self.batches_map[batch_name]['batched_data'].append(val)
            else:
                self.batches_map[batch_name]['batched_data'] = self.batches_map[batch_name]['batched_data'][-(self.batch_size - 1):] + [val]

        pass

    # TODO: handle case of not allowing to add a non-empty value!

    def _append_value_to_batch(self, batch_name, value_name, values_dict):
        val = values_dict[value_name]

        # Add the value (either FIFO style or at the end)
        if self.batch_size is None:
            self.batches_map[batch_name]['batched_data'].append(val)
        else:
            self.batches_map[batch_name]['batched_data'] = self.batches_map[batch_name]['batched_data'][-(self.batch_size - 1):] + [val]
        pass

    def _append_value_to_batch_array(self, batch_name, value_name, values_dict, start, end):
        val = values_dict[value_name]

        for i in range(start, end+1):
            if self.batch_size is None:
                self.batches_map[batch_name][i].append(val)
            else:
                self.batches_map[batch_name][i] = self.batches_map[batch_name][i][-(self.batch_size - 1):] + [val]
            pass

        pass

    def _is_value_empty(self, value):
        """
        If value is an array, return True/False depending on its length
        Otherwise, return True/False dependent on check if value is None
        """
        if type(value) == list:
            return len(value) == 0
        elif type(value) is None:
            return True

        return False
