from tkinter import *
import pandas as pd
import os
import csv
from tkinter import messagebox

# # attention, I used a semicolon as seperator
# Exp_list = pd.read_csv('U:\\Lab_2023\\Experiment_results\\QRAM\\20230717\\daily_experiment_comments.csv',
#                        usecols=lambda x: x != "Date")
# # iterate through all rows
# Valid_rows_str = []
# for row_index, row in Exp_list.iterrows():
#     # iterate through all elements in the row
#     if row['IgnoreValid'] == 'valid':
#         row_str = 'Time: %06d, ' % int(row['Time']) + ('with atoms, ' if row['Atoms'] else 'without atoms, ') + \
#                 '%s cycles, ' % row['Cycles'] + 'Comment: %s' % row['Comment']
#         Valid_rows_str.append([row_index, row_str])
#         print(row_str)

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.dirname(ROOT_DIR)
exp_path = []

Exp_tk = Tk()
Exp_tk.title('Experiment data load to dictionary')
# Exp_tk.geometry('750x400')


def _find_path_to_folder_containing_substring(substr, root_path):
    for dirname in os.listdir(root_path):
        if substr in dirname:
            return os.path.join(root_path, dirname)
    return None


def LoadList():
    # folder_path = os.path.join(DATA_PATH, exp_type.get(), exp_date.get())
    exp_data_path = os.path.join(DATA_PATH, exp_type.get())
    if not os.path.exists(exp_data_path):
        return messagebox.showwarning(title='Error!', message="Could not find a folder for experiment " + exp_type.get())

    exp_date_path = os.path.join(exp_data_path, exp_date.get())
    if not os.path.exists(exp_date_path):
        return messagebox.showwarning(title='Error!', message="Could not find in " + exp_type.get() +
                                                              " a folder for the date " + exp_date.get())
    filenames = next(os.walk(exp_date_path), (None, None, []))[2]
    Exp_list = pd.read_csv(os.path.join(exp_date_path, filenames[0]), usecols=lambda x: x != "Date")
    # iterate through all rows
    valid_rows_str = []
    for row_index, row in Exp_list.iterrows():
        # iterate through all elements in the row
        if row['IgnoreValid'] == 'valid':
            row_str = 'Time: %06d, ' % int(row['Time']) + ('with atoms, ' if row['Atoms'] else 'without atoms, ') + \
                      '%s cycles, ' % row['Cycles'] + 'Comment: %s' % row['Comment']
            valid_rows_str.append([row_index, row_str])
            print(row_str)

    for widget in frame2.winfo_children():
        widget.destroy()

    if len(valid_rows_str) > 0:
        for index, row in valid_rows_str:
            cb_exp.append(IntVar())
            cb_exp[-1].set(-1)
            Checkbutton(frame2, text=row, variable=cb_exp[-1], onvalue=index, offvalue=-1).pack(anchor='w')
    btn_load = Button(frame2, text='Load files', command=LoadFiles, padx=20, pady=5)
    btn_load.pack(side="bottom")

    return Exp_list


def isChecked():
    for cb in cb_exp:
        if cb.get() >= 0:
            Chosen_experiments.append(cb.get())


def LoadFiles():
    isChecked()
    if len(Chosen_experiments) > 0:
        str_files = ''
        for indx in Chosen_experiments:
            exp_path.append(_find_path_to_folder_containing_substring(str(Exp_list.Time[indx]), os.path.join(DATA_PATH, exp_type.get(), exp_date.get())))
        # data = Experiment_data_load.DictionaryBuilder(exp_type, exp_date, exp_time)
        # if data is None:
        #     exp_type, exp_date, exp_time = popupbox_inquiry()
        # else:
        #     Exp_dict = data.load_files_to_dict(data.exp_path)
        #     pymsgbox.alert(text='Experiment data is ready to use.', title='Success!')
        #     break
            str_files += f'The files from {exp_date.get()} at {Exp_list.Time[indx]} were loaded.\n'
            print(exp_type.get(), exp_date.get(), Exp_list.Time[indx])
        Exp_tk.destroy()
        return messagebox.showinfo('Load files', str_files)
    else:
        return messagebox.showwarning('Load files', 'Please Chose experiments to load!')


Chosen_experiments = []
cb_exp = []

frame1 = Label(Exp_tk, bg='#dddddd')
frame1.pack(side="top")
frame2 = LabelFrame(frame1, text='Valid Experiments List', padx=30, pady=10)
frame2.pack(side="bottom")

# Label(frame1, text='Experiment type').grid(row=0, column=0, padx=5, pady=5)
# Label(frame1, text='Experiment Date').grid(row=1, column=0, padx=5, pady=5)

Label(frame1, text='Experiment type').pack(side="left", fill="y")
exp_type = StringVar()
exp_type.set('QRAM')
exp_type_entry = Entry(frame1, textvariable=exp_type)
exp_type_entry.pack(side="left", fill="y")

Label(frame1, text='Experiment Date').pack(side="left", fill="y")
exp_date = StringVar()
exp_date.set('20230717')
exp_date_entry = Entry(frame1, textvariable=exp_date)
exp_date_entry.pack(side="left", fill="y")

btn_search = Button(frame1, text='Search', command=LoadList)
btn_search.pack(side="left")

Exp_list = LoadList()

Exp_tk.mainloop()

