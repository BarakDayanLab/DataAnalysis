import os
import glob
import pandas as pd
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import messagebox
from experiment_data import QRAMExperimentData
from exp_list import ExpListFrame
from config import DATA_PATH, EXP_TYPES

class ChooseExpFrame(tk.Frame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.master = master
        self.controller = controller
        self.exp_type = None
        self.exp_date = None
        self.search_results = []

        self.search_frame = SearchFrame(self, self.controller)
        self.search_frame.pack(side="top")

        self.exp_list_frame = ExpListFrame(self, self.controller, "multiple")
        self.exp_list_frame.pack(side="top", fill="both", expand=True)

        self.select_button = tk.Button(self, text="Select", command=self.select)
        self.select_button.pack(side="top")

    def update_search_results(self, exp_list):
        self.search_results = exp_list
        string_list = list(map(lambda exp: exp.desc, exp_list))
        self.exp_list_frame.update_exp_list(string_list)        

    def select(self):
        chosen_exps_idx = self.exp_list_frame.get_selection()
        exps = [self.search_results[idx] for idx in chosen_exps_idx]
        self.controller.set_chosen_exps(exps)


class SearchFrame(tk.Frame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.master = master
        self.controller = controller

        tk.Label(self, text="Experiment type:").grid(row=0, column=0)

        self.exp_type_entry = tk.Entry(self, 
                                       width=20,
                                       validate="focusout",
                                       validatecommand=self.validate_exp_type,
                                       highlightthickness=1)
        self.exp_type_entry.insert(tk.END, "QRAM")
        self.exp_type_entry.grid(row=0, column=1)

        self.exp_type_label = tk.Label(self, text="", foreground="red")
        self.exp_type_label.grid(row=1, column=1)

        tk.Label(self, text="Experiment date:").grid(row=0, column=2)
        self.exp_date_entry = tk.Entry(self, 
                                       width=20, 
                                       validate="focusout", 
                                       validatecommand=self.validate_exp_date,
                                       highlightthickness=1)
        self.exp_date_entry.insert(tk.END, "20230719")
        self.exp_date_entry.grid(row=0, column=3)

        self.exp_date_label = tk.Label(self, text="", foreground="red")
        self.exp_date_label.grid(row=1, column=3)

        self.search_button = tk.Button(self, text="Search", command=self.search)
        self.search_button.grid(row=0, column=4)

    def search(self):
        self.master.exp_type = self.exp_type_entry.get()
        self.master.exp_date = self.exp_date_entry.get()

        exp_date_path = os.path.join(DATA_PATH, self.master.exp_type, self.master.exp_date)
        
        exps_df = self.load_exps_csv(exp_date_path)
        exp_list = [QRAMExperimentData(self.master.exp_type, row) for row in exps_df.iterrows()]
        
        self.master.update_search_results(exp_list)

    def validate_exp_type(self):
        exp_type = self.exp_type_entry.get()
        exp_data_path = os.path.join(DATA_PATH, exp_type)
        if not os.path.exists(exp_data_path):
            self.exp_type_entry.config(highlightbackground="red")
            self.search_button.config(state="disabled")
            self.exp_type_label.config(text="Type does not exist!")
        else:
            self.exp_type_entry.config(highlightbackground="green")
            self.exp_type_label.config(text="")
            self.validate_exp_date()
        return True
    
    def validate_exp_date(self):
        exp_type = self.exp_type_entry.get()
        exp_date = self.exp_date_entry.get()
        exp_date_path = os.path.join(DATA_PATH, exp_type, exp_date)
        if not os.path.exists(exp_date_path):
            self.exp_date_entry.config(highlightbackground="red")
            self.search_button.config(state="disabled")
            self.exp_date_label.config(text="Date does not exist!")
        else:
            self.search_button.config(state="active")
            self.exp_date_entry.config(highlightbackground="green")
            self.exp_date_label.config(text="")
        return True

        
    @staticmethod
    def load_exps_csv(exp_date_path):
        csv_filepath = glob.glob(os.path.join(exp_date_path, "*.csv"))
        if len(csv_filepath) != 1:
            messagebox.showerror(title="Error!", 
                                message=f"Could not find the csv file in \n{exp_date_path}")
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
    

class NewSearchFrame(tk.Frame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.master = master
        self.controller = controller

        tk.Label(self, text="Experiment type:").grid(row=0, column=0)

        exp_types = [d for d in os.listdir(DATA_PATH) if os.path.isdir(os.path.join(DATA_PATH)) and d in EXP_TYPES] 
        self.exp_type = tk.StringVar(value=exp_types[0])
        self.exp_type_comb = ttk.Combobox(self, values=exp_types, textvariable=self.exp_type).grid(row=0, column=1)

        tk.Label(self, text="Experiment date:").grid(row=0, column=2)
        self.exp_date_entry = tk.Entry(self, 
                                       width=20, 
                                       validate="focusout")
        
        exp_dates = tk.StringVar([])
        exp_dates.set([d for d in os.listdir(DATA_PATH) if os.path.isdir(os.path.join(DATA_PATH)) and d in EXP_TYPES])
        exp_dates = [d for d in os.listdir(DATA_PATH) if os.path.isdir(os.path.join(DATA_PATH)) and d in EXP_TYPES] 
        self.exp_date_entry.insert(tk.END, "20230719")
        self.exp_date_entry.grid(row=0, column=3)

        self.exp_date_label = tk.Label(self, text="", foreground="red")
        self.exp_date_label.grid(row=1, column=3)

        self.search_button = tk.Button(self, text="Search", width=7, command=self.search)
        self.search_button.grid(row=0, column=4)

    
    def search(self):
        exp_type = self.exp_type.get()#self.exp_type_entry.get()
        exp_date = self.exp_date_entry.get()
        print(exp_type, exp_date)

        exp_date_path = os.path.join(DATA_PATH, exp_type, exp_date)
        
        exps_df = self.load_exps_csv(exp_date_path)
        exp_list = self.parse_df_to_list(exps_df)
        # self.
        self.controller.exp_list_frame.update_exp_list(exp_list)

        
    @staticmethod
    def load_exps_csv(exp_date_path):
        csv_filepath = glob.glob(os.path.join(exp_date_path, "*.csv"))
        if len(csv_filepath) != 1:
            messagebox.showerror(title="Error!", 
                                message=f"Could not find the csv file in \n{exp_date_path}")
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
        
