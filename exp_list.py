import tkinter as tk

class ExpListFrame(tk.Frame):
    def __init__(self, master, controller, selectmode, selection_callback=None):
        super().__init__(master)
        self.master = master
        self.controller = controller
        self.selectmode = selectmode

        self.exp_list = tk.StringVar()
        self.exp_list.set([])

        tk.Label(self, text="Experiment list:").pack(side="top")
        self.exp_listbox = tk.Listbox(self, listvariable=self.exp_list, selectmode=selectmode)
        if selection_callback is not None:
            self.exp_listbox.bind("<<ListboxSelect>>", selection_callback)
        self.exp_listbox.pack(side="top", fill="both", expand=True)

    def update_exp_list(self, new_exp_list):
        self.exp_list.set(new_exp_list)

    def get_selection(self):
        return list(self.exp_listbox.curselection())