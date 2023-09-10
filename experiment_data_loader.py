import os
from typing import Any
import numpy as np
from glob import glob
from tkinter import messagebox

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.dirname(ROOT_DIR)


# class MultipleExpDataLoader:
#     def __init__(self):
#         self.exp_type = None
#         self.exps_df = None
#         self.exps_data = None 
    
#     @staticmethod
#     def load_exp_data(exp_path):
#         data_dir = os.path.join(DATA_PATH, exp_path)
#         npz_files = glob(os.path.join(data_dir, "**", "*.npz"), recursive=True)
#         if len(npz_files) == 0:
#             messagebox.showerror(title='Error!', 
#                                 message=f"Could not find in npz files in {exp_path}")
#             return
        
#         exps_dict = {}
#         for npz_file in npz_files:
#             exp_dict = np.load(npz_file, allow_pickle=True)["arr_0"]

#             keys = npz_file[len(data_dir+os.sep):-4].split(os.sep)
#             print(keys)
#             for key in keys[::-1]:
#                 exp_dict = {key: exp_dict}
            
#             exps_dict = {**exps_dict, **exp_dict}

#         return exps_dict   

#     def load_exps(self, exp_type, exps_df):
#         self.exp_type = exp_type
#         self.exps_df = exps_df

#         exp_date_list = self.exps_df["Date"].tolist()
#         exp_time_list = self.exps_df["Time"].tolist()
#         exp_paths = [os.path.join(self.exp_type, str(date), f"{time}_Photon_TimeTags") for date, time in zip(exp_date_list, exp_time_list)]
#         self.exps_data = list(map(self.load_exp_data, exp_paths))

class MultipleExpDataLoader:
    def __init__(self, exp_type, exps_df):
        self.exp_type = exp_type
        self.exps_df = exps_df

        df_list = self.exps_df.apply(lambda x: (str(x["Date"]), str(x["Time"])), axis=1).tolist()
        self.exps_data: list[ExpDataLoader] = list(map(self.df_to_exp_list, df_list))

    def df_to_exp_list(self, val):
        date, time = val
        return ExpDataLoader(self.exp_type, date, time)

class ExpDataLoader:
    def __init__(self, exp_type, exp_date, exp_time):
        self.exp_type = exp_type
        self.exp_date = exp_date
        self.exp_time = exp_time
        self.exp_path = os.path.join(DATA_PATH, self.exp_type, self.exp_date, f"{self.exp_time}_Photon_TimeTags")
        self.curr_attr_path = self.exp_path


    def __getattr__(self, name):
        self.curr_attr_path = os.path.join(self.curr_attr_path, name)
        # print(self.curr_attr_path, name)

        if os.path.isdir(self.curr_attr_path):
            return self
        else:
            file_data = self._load_npz_file(f"{self.curr_attr_path}.npz")
            self.curr_attr_path = self.exp_path
            return file_data


    def _load_npz_file(self, npz_file):
        return np.load(npz_file, allow_pickle=True)["arr_0"]


if __name__ == '__main__':
    print(ExpDataLoader("QRAM", "20230719", "010507").BalancingRes.MZ_balancing_Phase_Correction_min_vec_batch.shape)