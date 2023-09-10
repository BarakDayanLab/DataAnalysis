import os
import glob
import pandas as pd
import tkinter as tk
from tkinter import messagebox
import tkinter.ttk as ttk
from experiment_data_loader import MultipleExpDataLoader
from filter_echos import plot_data

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.dirname(ROOT_DIR)

class ChooseExpGui(tk.Tk):
    def __init__(self):
        super().__init__()
        self.data_loader = None

        self.wm_title("Experiment data load to dictionary")
        self.main_container = tk.Frame(self)
        self.main_container.pack(side="top", fill="both", expand=True)

        self.search_frame = SearchFrame(self.main_container, self)
        self.search_frame.pack(fill="both", expand=True)

        self.exp_list_frame = ExpListFrame(self.main_container, self)
        self.exp_list_frame.pack(side="top", fill="both", expand=True)

class SearchFrame(tk.Frame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.master = master
        self.controller = controller
        self.exp_type = None
        self.exp_date = None
        self.exp_df = None

        tk.Label(self, text="Experiment type:").pack(side="left", expand=True)

        self.exp_type_entry = tk.Entry(self, width=20)
        self.exp_type_entry.insert(tk.END, "QRAM")
        self.exp_type_entry.pack(side="left", expand=True)

        tk.Label(self, text="Experiment date:").pack(side="left", expand=True)
        self.exp_date_entry = tk.Entry(self, width=20)
        self.exp_date_entry.insert(tk.END, "20230719")
        self.exp_date_entry.pack(side="left", expand=True)

        self.search_button = tk.Button(self, text="Search", width=7, command=self.search)
        self.search_button.pack(side="left", expand=True)

    def search(self):
        exp_type = self.exp_type_entry.get()
        exp_date = self.exp_date_entry.get()

        exp_data_path = os.path.join(DATA_PATH, exp_type)
        if not os.path.exists(exp_data_path):
            self.controller.exp_list_frame.update_exp_list([])
            self.exp_df = None
            messagebox.showerror(title="Error!", 
                                message=f"Could not find a folder for experiment {exp_type}")
            return 
        
        exp_date_path = os.path.join(exp_data_path, exp_date)
        if not os.path.exists(exp_date_path):
            self.controller.exp_list_frame.update_exp_list([])
            self.exp_df = None
            messagebox.showerror(title='Error!', 
                                message=f"Could not find in {exp_type} a folder for the date {exp_date}")
            return 
        
        self.exp_df = self.load_exps_csv(exp_date_path)
        exp_list = self.parse_df_to_list(self.exp_df)
        self.controller.exp_list_frame.update_exp_list(exp_list)
        
    @staticmethod
    def load_exps_csv(exp_date_path):
        csv_filepath = glob.glob(os.path.join(exp_date_path, "*.csv"))
        if len(csv_filepath) != 1:
            messagebox.showerror(title="Error!", 
                                message=f"Could not find the csv file in {exp_date_path}")
            return
        
        csv_filepath = csv_filepath[0]

        exp_df = pd.read_csv(csv_filepath)
        valid_exp_df = exp_df.query("IgnoreValid == 'valid'")
        return valid_exp_df
    
    @staticmethod
    def parse_df_to_list(df):
        parse_row = lambda row: f"Time: {row['Time']}, {'with atoms' if row['Atoms'] else 'without atoms'}, {row['Cycles']} cycles, Comment: {row['Comment']}"
        exp_list = df.apply(parse_row, axis=1).tolist()
        return exp_list


class ExpListFrame(tk.Frame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.master = master
        self.controller = controller

        self.exp_list = tk.StringVar()
        self.exp_list.set([])

        tk.Label(self, text="Experiment list:").pack(side="top", expand=True)
        self.exp_listbox = tk.Listbox(self, listvariable=self.exp_list, width=50, height=20, selectmode="multiple")
        self.exp_listbox.pack(side="top", fill="both", expand=True)

        self.select_button = tk.Button(self, text="Select", width=7, command=self.select)
        self.select_button.pack(side="top")

    def update_exp_list(self, new_exp_list):
        self.exp_list.set(new_exp_list)

    def select(self):
        chosen_exps_idx = list(self.exp_listbox.curselection())
        chosen_exps = self.controller.search_frame.exp_df.iloc[chosen_exps_idx]
        exp_type = self.controller.search_frame.exp_type_entry.get()
        self.controller.data_loader = MultipleExpDataLoader(exp_type, chosen_exps)
        plot_data(self.controller.data_loader.exps_data[0].North_timetags_folded_to_seq)
        # self.controller.data_loader = MultipleExpDataLoader()
        # self.controller.data_loader.load_exps(exp_type, chosen_exps)
        # print(self.controller.data_loader.exps_data[0].keys())


    



if __name__ == "__main__":
    root = ChooseExpGui()
    root.mainloop()