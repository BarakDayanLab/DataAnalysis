import os
import csv
import json
import time
import errno
import numpy as np
import shutil
from scipy.io import savemat
import collections.abc
import matplotlib.pyplot as plt
from Utilities.Utils import Utils


class BDResults:

    def __init__(self, json_map_path=None, version=None, logger=None):

        self.strict = True

        if logger is not None:
            self.logger = logger

        # Load results map from json file
        try:
            the_path = '.' if json_map_path is None else json_map_path
            the_file = os.path.join(the_path, 'results_map.json')
            f = open(the_file)
            self.results_map = json.load(f)
            f.close()
        except Exception as err:
            self.logger.error(f'Unable to open results_map.json. Reason: {err}')
            raise Exception(f'Unable to open results_map.json. Reason: {err}')

        self.strict = self.results_map['strict'] if 'strict' in self.results_map else None

        # Resolve root folder and create it
        if 'root' not in self.results_map:
            self.logger.error("Cannot find 'root' in results_map.json.")
            raise Exception("Cannot find 'root' in results_map.json.")

        if 'files' not in self.results_map:
            self.logger.error("Cannot find 'files' in results_map.json.")
            raise Exception("Cannot find 'files' in results_map.json.")

        # We initialize the folders map with None, since it is not relevant until someone
        # invokes the create_folders() - so they are actually resolved and created.
        self.folders = None

        self.experiment_run_folder = None

        # Perform version check
        if version is not None:
            if 'version' not in self.results_map:
                self.logger.error('results map file has no version mentioned')
            elif self.results_map['version'] != version:
                self.logger.error(f"results map version mismatch! (is {self.results_map['version']} but should be {version})")
                raise Exception(f"results map version mismatch! (is {self.results_map['version']} but should be {version})")
        pass

    def _create_folder(self, path, allow_existing=True):
        # Create a folder and check if there's already one there
        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
        pass

    def _save_file(self, file_entity, path, data_pool):
        data_attr_to_save = file_entity['data'] if 'data' in file_entity else None
        data_to_save = data_pool[data_attr_to_save] if data_attr_to_save in data_pool else None
        name_to_save = file_entity['name'] if 'name' in file_entity else None
        should_append = file_entity['append'] if 'append' in file_entity else None
        header = file_entity['header'] if 'header' in file_entity else None
        types = file_entity['type'] if 'type' in file_entity else None

        if not self._is_array(types):
            types = [types]
        for _type in types:
            if _type == 'csv':
                self._save_csv_file(data_to_save, name_to_save, path, should_append, header)
            elif _type == 'csv_writer':
                self._save_csv_writer_file(data_to_save, path)
            elif _type == 'txt':
                self._save_txt_file(data_to_save, path, should_append)
            elif _type == 'mat':
                self._save_matlab_file(data_to_save, name_to_save, path)
            elif _type == 'npz':
                self._save_numpy_file(data_to_save, path)
            elif _type == 'plt':
                self._save_figure(path)
            elif _type == 'json':
                self._save_json(data_to_save, path)
            else:
                self.logger.error(f'Invalid file type- {_type}. Check your json configuration.')

    def _save_figure(self, path):
        # Ensure there is the ".png" extension
        if not path.endswith('.png'):
            path = path + '.png'
        plt.savefig(path, bbox_inches='tight')

    def _save_txt_file(self, data, path, should_append):
        # Ensure there is the ".txt" extension
        if not path.endswith('.txt'):
            path = path + '.txt'
        mode = "a" if should_append else "w"
        with open(path, mode) as text_file:
            text_file.write(data + '\n')
        pass

    # TODO: we may want to change this function to write data objects without the need to format it before
    def _save_csv_file(self, data, name, path, should_append, header):
        # Ensure there is the ".csv" extension
        if not path.endswith('.csv'):
            path = path + '.csv'
        mode = "a" if should_append else "w"

        # If this is the first time, we need to write the header
        should_write_header = header and not os.path.exists(path)

        with open(path, mode) as text_file:
            if should_write_header:
                text_file.write(header + '\n')
            text_file.write(data + '\n')
        pass

    def _save_csv_writer_file(self, data, path):
        """
        This function assumes the data object is of this format:
        data = {
            "header_1": [ val1, val2, val3, ...]
            "header_2": [ val1, val2, val3, ...]
            ...
            "header_N": [ val1, val2, val3, ...]
        }
        """
        # Ensure there is the ".csv" extension
        if not path.endswith('.csv'):
            path = path + '.csv'

        try:
            with open(path, "w", newline='', encoding='utf-8') as outfile:
                writer = csv.writer(outfile)
                # Write the headers
                writer.writerow(data.keys())
                # Write the data
                writer.writerows(zip(*data.values()))
        except IOError as e:
            self.logger.error(f'I/O error in saveLinesAsCSV {e}')

        pass
    def _save_numpy_file(self, data, path):
        # Ensure there is the ".npz" extension
        if not path.endswith('.npz'):
            path = path + '.npz'
        # Save the file
        np.savez(path, data)
        pass

    def _save_matlab_file(self, data, name, path):
        # Ensure there is the ".mat" extension
        if not path.endswith('.mat'):
            path = path + '.mat'
        # Save the file
        savemat(path, mdict={name: data})
        pass

    def _save_json(self, data, path):
        # Ensure there is the ".json" extension
        if not path.endswith('.json'):
            path = path + '.json'
        # Save the file
        with open(path, 'w') as file:
            json.dump(data, file, indent=4)
        pass

    def _is_array(self, obj):
        """
        Checks if a given object is an array or np.array (not a scalar, float or string...)
        :param obj: Object to check
        :return: True/False
        """
        return not isinstance(obj, str) and isinstance(obj, (collections.abc.Sequence, np.ndarray))

    def _resolve_env(self, root_name):
        """
        Given a specific root name, this method first looks for such a root with the current ENV
        If not found, it fetches the root without the ENV
        """
        env = self.results_map['env']
        root_name_with_env = f'{root_name}_{env}'
        if root_name_with_env in self.results_map:
            return self.results_map[root_name_with_env]

        if root_name not in self.results_map:
            self.logger.error(f'Missing root {root_name} in results_map.json. Please check.')

        return self.results_map[root_name]


    def _resolve_parameterized(self, root_pattern):
        current_date_time = time.strftime("%Y%m%d_%H%M%S")
        current_date = time.strftime("%Y%m%d")
        current_time = time.strftime("%H%M%S")

        root_pattern = root_pattern.replace('{current_date_time}', current_date_time)
        root_pattern = root_pattern.replace('{current_date}', current_date)
        root_pattern = root_pattern.replace('{current_time}', current_time)

        return root_pattern

    def get_root(self):
        resolved_env = self._resolve_env('root')
        resolved_path = self._resolve_parameterized(resolved_env)
        return resolved_path

    def get_experiment_root(self):
        resolved_env = self._resolve_env('experiment_root')
        resolved_path = self._resolve_parameterized(resolved_env)
        return resolved_path

    def get_all_experiments_root(self):
        resolved_env = self._resolve_env('all_experiments_root')
        resolved_path = self._resolve_parameterized(resolved_env)
        return resolved_path

    def get_custom_root(self, custom_root):
        resolved_env = self._resolve_env(custom_root)
        resolved_path = self._resolve_parameterized(resolved_env)
        return resolved_path

    def get_folder_path(self, name):
        if self.folders is None:
            raise Exception('You must invoke create_folders() first to resolve/create the relevant folders')

        return self.folders[name]

    def get_experiment_run_folder(self):
        return self.experiment_run_folder

    def get_sequence_folder(self, sequence_definitions):
        """
        Return the name of the folder according to the iteration and sequence-step
        Example: 20240221_140500_Iteration 1_Sequence 1__No Atoms
        """
        experiment_run_folder = self.get_experiment_run_folder()

        iteration = sequence_definitions['current_iteration']
        sequence_step = sequence_definitions['sequence_step']
        name = sequence_definitions['sequence'][sequence_step]['name']

        sequence_folder_path = f'Iter_{iteration+1}_Seq_{sequence_step+1}__{name}'

        # Construct base and the sequence
        sequence_folder_path = os.path.join(experiment_run_folder, sequence_folder_path)

        return sequence_folder_path

    def create_experiment_run_folder(self):
        """
        Resolve the SPECIFIC experiment RUN folder and create it (if does not exist)
        This is the folder with the time-stamp
        It uses the "root" entry in the results_map.json
        """
        self.experiment_run_folder = self._resolve_parameterized(self.results_map['root'])

        if not os.path.exists(self.experiment_run_folder):
            os.makedirs(self.experiment_run_folder, exist_ok=True)

        return self.experiment_run_folder

    def save_results(self, data_pool, resolved_path):

        # Resolve the file names and folders in results map
        for entity in self.results_map['files']:
            if 'file_name' in entity:
                entity['file_name'] = self._resolve_parameterized(entity['file_name'])
            if 'folder' in entity:
                entity['folder'] = self._resolve_parameterized(entity['folder'])

        # Iterate over all result entities we need to save
        for entity in self.results_map['files']:
            entity_keys = list(entity.keys())
            # Is this a comment line?
            if len(entity_keys) == 1 and entity_keys[0].startswith('#'):
                continue

            if 'skip' in entity and entity['skip']:
                continue

            # Sanity check on type
            if 'type' not in entity:
                self.logger.error("Missing 'type' entity in 'file' attribute")
                continue

            # Is this a data vector we need to save (not a plot)
            if entity['type'] != 'plt':
                if 'data' not in entity:
                    self.logger.error("Missing 'data' entity in 'file' attribute")
                    continue
                if entity['data'] not in data_pool:
                    self.logger.error(f"Cannot find {entity['data']} in the results data")
                    continue

            # Resolve the root folder - with current time/date
            path = resolved_path
            if 'folder' in entity:
                path = os.path.join(resolved_path, entity['folder'])

            # Create folder if it does not exist
            if not os.path.exists(path):
                os.makedirs(path, exist_ok=True)

            # Add the file name
            file_path = os.path.join(path, entity['file_name'])

            start_save_time = time.time()

            # Save the file
            try:
                self._save_file(entity, file_path, data_pool)
            except Exception as err:
                self.logger.error(f'Failed to save file {file_path}. {err}')

            self.logger.info(f'Saving {file_path} - Done! (took {time.time() - start_save_time} sec)')

        return resolved_path

    def create_folders(self):
        """
        This method will create all folders defined in results_map.json
        NOTE - it will do the path resolving and folder creation ONLY when it's called.
               So, if there are {current_date} and/or {current_time} fields, they will get the NOW time
               Until this method is not called, one cannot call get_folder_path(...) - they are simply not resolve/created
        """
        resolved_path = self._resolve_parameterized(self.results_map['root'])

        self.folders = {}

        # Resolve and insert to dictionary the various root folders first
        for root_name in ['root', 'experiment_root', 'all_experiments_root']:
            if root_name in self.results_map:
                self.folders[root_name] = self._resolve_parameterized(self.results_map[root_name])

        # Resolve the file names and folders in results map
        for entity in self.results_map['folders']:
            if 'name' in entity:
                entity['name'] = self._resolve_parameterized(entity['name'])

            # Create folder if it does not exist
            path = os.path.join(resolved_path, entity['name'])
            if not os.path.exists(path):
                os.makedirs(path, exist_ok=True)

            # Create a resolved entry in the folders map
            self.folders[entity['name']] = path

        pass

    def copy_files(self, source, destination, opt_in_filter=None, opt_out_filter=None, create_folder=False):
        """
        Copy files from source to destination
        - source - source folder of files
        - destination - destination folder to copy files
        - opt_in_filter - Only include files that contain this filter
        - opt_out_filter - Exclude files that contain this filter
        - create_folder - if source/destination folders do not exist, create them
        """

        # Create folders if they do not exist
        if create_folder:
            if not os.path.exists(source):
                os.makedirs(source)
            if not os.path.exists(destination):
                os.makedirs(destination)

        try:
            files_to_copy = Utils.get_files_in_path(source, opt_in_filter=opt_in_filter, opt_out_filter=opt_out_filter, return_full_path=False)
            for file in files_to_copy:
                src = os.path.join(source, file)
                dst = os.path.join(destination, file)
                shutil.copy2(src, dst)
        except Exception as err:
            self.logger.error(f'Failed to copy files from {source} to {destination}. {err}')

        pass

    def copy_folder(self, source, destination, create_destination_folder=False):
        """
        Copy recursively entire folder from source to destination. Resolve destination before copy.
        If source is not a folder, it ignores it (does not copy!)
        """

        try:
            source = self._resolve_parameterized(source)
            destination = self._resolve_parameterized(destination)

            # If destination folder does not exist, create it:
            if create_destination_folder and not os.path.exists(destination):
                os.makedirs(destination, exist_ok=True)

            destination_2 = shutil.copytree(source, destination)
        except OSError as err:
            # error caused if the source was not a directory
            if err.errno == errno.ENOTDIR:
                self.logger.error(f'Destination {destination} is not a folder. {err}')
            else:
                self.logger.error(f'Copy Folder failed (src: {source}, dst: {destination}). {err}')

    def is_network_drive_available(self, network_path):
        success = True

        try:
            network_file = os.path.join(network_path, '__test_file__.txt')
            # Try to write a file
            with open(network_file, "w") as file:
                file.write("Test File - checking access to network drive")
            # Remove the file
            os.remove(network_file)
        except Exception as err:
            success = False

        return success

    @staticmethod
    def unit_tests():

        # Open a figure, so we can test saving a plot
        plt.figure()
        plt.plot(np.array([1, 2, 3, 4]), label='Something 1')
        plt.legend()
        plt.title('This is a title')
        plt.show(block=False)
        plt.pause(0.01)

        # Create the data we want to save
        data = {
            "comments": {
                "json_key_1": "key_1",
                "json_key_2": "key_2",
                "array": [{"a": 1}, {"b": 2}, {"c": 3}],
                "flag": True
            },
            "exp_comment": "transit condition: minimum time between reflection = 2000with at least 3 reflections ; reflection threshold: 0.2MCounts / sec",
            "exp_comment_csv": {
                "date": "",
                "time": "",
                "success": "",
                "with_atoms": "",
                "counter": "",
                "comment": ""
            },
            "daily_experiment_comments": "20230831,102921,ignore,True,2,ignore",
            "input_vector": np.array([1, 2, 3, 4, 5]),
            "output_vector": [6, 7, 8, 9],
            "test_data_for_csv_writer": {
                "header_a": [1, 2, 3],
                "header_b": [4, 5, 6],
                "header_z": [-1, -3, -9]
            }
        }

        # Initiate and test
        bd_results = BDResults(version="0.1")
        bd_results.save_results(data)

        pass

if __name__ == "__main__":
    BDResults.unit_tests()

