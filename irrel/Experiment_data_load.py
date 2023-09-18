import os
import shutil
from os import listdir
from os.path import isfile, join
from os import walk
import pandas as pd
import numpy as np
import pymsgbox


class DictionaryBuilder:

    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_PATH = os.path.dirname(ROOT_DIR)

    def first_index_containing_substring(self, the_list, substring):
        for i, s in enumerate(the_list):
            if substring in s:
                return i
        return -1

    def _find_background(self, indx):
        if (indx > 0) and not self.daily_experiment_comments.Atoms[indx - 1] and \
                (self.daily_experiment_comments.IgnoreValid[indx - 1] == 'valid') and \
                (abs(int(self.daily_experiment_comments.Time[indx]) - int(self.daily_experiment_comments.Time[indx])) //
                 1e4 < 2):
            self.runs_data_and_background_duas.append([self.daily_experiment_comments.Time[indx],
                                                       self.daily_experiment_comments.Time[indx - 1]])
        elif not self.daily_experiment_comments.Atoms[indx + 1] and \
                (self.daily_experiment_comments.IgnoreValid[indx + 1] == 'valid') and \
                (abs(int(self.daily_experiment_comments.Time[indx]) - int(self.daily_experiment_comments.Time[indx])) //
                 1e4 < 2):
            self.runs_data_and_background_duas.append([str(self.daily_experiment_comments.Time[indx]),
                                                       str(self.daily_experiment_comments.Time[indx + 1])])

    def _sorting_dirs_to_with_and_without_atom(self):
        folder_list = os.listdir(self.exp_date_path)
        filenames = next(walk(self.exp_date_path), (None, None, []))[2]
        self.daily_experiment_comments = pd.read_csv(os.path.join(self.exp_date_path, filenames[0]))
        self.dirname_with_atoms = os.path.join(self.exp_date_path, 'with_Atoms')
        self.dirname_without_atoms = os.path.join(self.exp_date_path, 'without_Atoms')
        if not os.path.exists(self.dirname_with_atoms):
            os.makedirs(self.dirname_with_atoms)
        if not os.path.exists(self.dirname_without_atoms):
            os.makedirs(self.dirname_without_atoms)

        for indx, time in enumerate(self.daily_experiment_comments.Time):
            folder_indx = self.first_index_containing_substring(folder_list, str(time))
            if folder_indx >= 0 and self.daily_experiment_comments.IgnoreValid[indx] == 'valid':
                if self.daily_experiment_comments.Atoms[indx]:
                    shutil.move(os.path.join(self.exp_date_path, folder_list[folder_indx]),
                                os.path.join(self.dirname_with_atoms, folder_list[folder_indx]))
                    self._find_background(indx)
                else:
                    shutil.move(os.path.join(self.exp_date_path, folder_list[folder_indx]),
                                os.path.join(self.dirname_without_atoms, folder_list[folder_indx]))

        if len(self.runs_data_and_background_duas) > 0:
            np.savez(os.path.join(self.exp_date_path, 'runs_data_and_background_duas.npz'), self.runs_data_and_background_duas)

    def _build_dict(self):
        self.Experiments_data = {
            self.exp_date: {}
        }
        for subdir in next(walk(self.exp_date_path))[1]:
            if subdir.split('_')[1] == 'Atoms':
                self.Experiments_data[self.exp_date][subdir] = self._load_files_to_dict(
                    os.path.join(self.exp_date_path, subdir))

    def _build_dict_dua(self, run_data_and_background_dua):
        self.Experiments_data = {
            self.exp_date: {
            }
        }
        for subdir in next(walk(os.path.join(self.exp_date_path, 'with_Atoms')))[1]:
            if run_data_and_background_dua[0] in subdir:
                self.Experiments_data[self.exp_date]['with_atoms'] = self._load_files_to_dict_dua(
                    os.path.join(self.exp_date_path, 'with_Atoms', subdir))
        for subdir in next(walk(os.path.join(self.exp_date_path, 'without_Atoms')))[1]:
            if run_data_and_background_dua[1] in subdir:
                self.Experiments_data[self.exp_date]['without_atoms'] = self._load_files_to_dict_dua(
                    os.path.join(self.exp_date_path, 'without_Atoms', subdir))

    def _load_files_to_dict_dua(self, folder_path):
        folder = {}
        itr = next(walk(folder_path))
        if len(itr[1]) > 0:
            for dirname in itr[1]:
                try:
                    folder[dirname] = self._load_files_to_dict_dua(os.path.join(folder_path, dirname))
                except Exception as err:
                    pass
        # else:
        for filename in itr[2]:
            if filename.split('.')[1] == 'npz':
                try:
                    folder[filename.split('.')[0]] = \
                        np.load(os.path.join(folder_path, filename), allow_pickle=True)["arr_0"]
                except Exception as err:
                    pass
        return folder

    def _find_path_to_folder_containing_substring(self, substr, root_path):
        for dirname in os.listdir(root_path):
            if substr in dirname:
                return os.path.join(root_path, dirname)
        return None

    def load_files_to_dict(self, folder_path):
        folder = {}
        itr = next(walk(folder_path))
        if len(itr[1]) > 0:
            for dirname in itr[1]:
                try:
                    folder[dirname] = self.load_files_to_dict(os.path.join(folder_path, dirname))
                except Exception as err:
                    pass
        # else:
        for filename in itr[2]:
            if filename.split('.')[1] == 'npz':
                try:
                    folder[filename.split('.')[0]] = \
                        np.load(os.path.join(folder_path, filename), allow_pickle=True)["arr_0"]
                except Exception as err:
                    pass
        return folder

    def _load_files(self):
        pass


    # Class's constructor
    def __init__(self, exp_type='QRAM', exp_date='20230719', exp_time=None):

        self.exp_data_path = os.path.join(self.DATA_PATH, exp_type)
        if not os.path.exists(self.exp_data_path):
            pymsgbox.alert(text="Could not find a folder for experiment " + exp_type, title='Error!')
            return None
        
        self.exp_date = exp_date
        self.exp_date_path = os.path.join(self.exp_data_path, self.exp_date)
        if not os.path.exists(self.exp_date_path):
            pymsgbox.alert(text="Could not find in " + exp_type + "a folder for the date " + exp_date, title='Error!')
            return None

        if exp_time is None:
            self.exp_path = self.exp_date_path
        else:
            self.exp_time = exp_time
            self.exp_path = self._find_path_to_folder_containing_substring(self.exp_time, self.exp_date_path)
        if self.exp_path is None:
            pymsgbox.alert(text="Could not find for " + exp_type + " experiment at " + exp_date +
                                " a folder for the requested time - " + exp_time, title='Error!')
            return None

        # self.Exp_dict = self._load_files_to_dict(self.exp_path)
        # self.runs_data_and_background_duas = []

if __name__ == '__main__':

        dict = DictionaryBuilder()
        print(dict.DATA_PATH)
