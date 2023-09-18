import os
import numpy as np
from glob import glob
from collections import UserDict
import utilities as utils

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.dirname(ROOT_DIR)

"""
TODO: add errors and warnings
"""

class LazyDict(UserDict):
    def __init__(self, dict, base_path):
        super().__init__(dict)
        self.base_path = base_path
        self.current_file_path = ""

    def __getitem__(self, key):
        if self.current_file_path == "":
            val = self.data.get(key)
            if val is not None:
                return val
         
        self.current_file_path = os.path.join(self.current_file_path, key)
        if os.path.isdir(os.path.join(self.base_path, self.current_file_path)):
            return self
        
        attrs = self.current_file_path.split(os.sep)
        file_data = self._load_npz_file(os.path.join(self.base_path, f"{self.current_file_path}.npz"))
        for idx, attr in enumerate(attrs[1:][::-1]):
            dict_base_path = os.path.join(self.base_path, os.sep.join(attrs[:len(attrs)-idx-1]))
            file_data = LazyDict({attr: file_data}, dict_base_path)
        self.data[attrs[0]] = file_data
        self.current_file_path = ""
        return file_data

    def _load_npz_file(self, npz_file):
        return np.load(npz_file, allow_pickle=True)["arr_0"]
    
    

class QRAMExperimentData:
    def __init__(self, exp_type, exp_date, exp_time, valid, atoms, cycles, comment):
        self.type = exp_type
        self.date = exp_date
        self.time = exp_time
        self.valid = valid
        self.atoms = atoms
        self.cycles = cycles
        self.comment = comment

        self.path = os.path.join(DATA_PATH, self.type, self.date, f"{self.time}_Photon_TimeTags")
        self.files = LazyDict({}, self.path)

    def __init__(self, exp_type, df_row):
        self.type = exp_type
        self.exp_idx, row = df_row
        self.date = str(row["Date"])
        self.time = str(row["Time"])
        self.valid = row["IgnoreValid"]
        self.atoms = row["Atoms"]
        self.cycles = row["Cycles"]
        self.comment = row["Comment"]

        self.path = os.path.join(DATA_PATH, self.type, self.date, f"{self.time}_Photon_TimeTags")
        self.files = LazyDict({}, self.path)

    @property
    def title(self):
        return f"{self.parsed_date} at {self.parsed_time}"

    @property
    def parsed_date(self):
        """
        parse date from the format yyyymmdd to dd/mm/yyyy
        """
        return f"{self.date[6:8]}/{self.date[4:6]}/{self.date[:4]}"
    
    @property
    def parsed_time(self):
        """
        parse time from the format hhmmss to hh:mm:ss
        """
        return f"{self.time[:2]}:{self.time[2:4]}:{self.time[4:6]}"
    
    @property
    def desc(self):
        return f"{self.parsed_date} at {self.parsed_time}, {'with atoms' if self.atoms else 'without atoms'}, {self.cycles} cycles, Comment: {self.comment}"

if __name__ == '__main__':
    dicti = LazyDict({}, os.path.join(DATA_PATH, "QRAM", "20230719", "010507_Photon_TimeTags"))

    print("")