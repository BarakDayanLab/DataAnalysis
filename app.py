import tkinter as tk
from search_exp import ChooseExpFrame
from view_exp import SelectedExpsFrame


class MainWindow(tk.Tk):
    def __init__(self):
        super().__init__()
        self.qram_exps = []

        self.geometry("1200x800")

        self.wm_title("Experiment data load to dictionary")
        self.main_container = tk.Frame(self)
        self.main_container.pack(side="top", fill="both", expand=True)

        self.choose_exp_frame = ChooseExpFrame(self.main_container, self)
        self.choose_exp_frame.pack(side="top", fill="both", expand=True)

        self.selected_exp_frame = SelectedExpsFrame(self.main_container, self)

    def set_chosen_exps(self, exps):
        self.qram_exps = exps
        self.show_selected_exp()  

    def show_selected_exp(self):
        self.choose_exp_frame.pack_forget()
        self.selected_exp_frame.pack(side="top", fill="both", expand=True)
        self.selected_exp_frame.update_exp_list()

    def back_to_search(self):
        self.selected_exp_frame.pack_forget()
        self.choose_exp_frame.pack(side="top", fill="both", expand=True)

if __name__ == "__main__":
    app = MainWindow()
    app.mainloop()