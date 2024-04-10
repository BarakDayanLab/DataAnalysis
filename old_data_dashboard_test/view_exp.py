import tkinter as tk
from exp_list import ExpListFrame
from choose_plot_popup import ChoosePlotPopup
from data_visualization import PlotManagerFrame, QRAMDataVisualization


class SelectedExpsFrame(tk.Frame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.master = master
        self.controller = controller
        self.open_popup = None

        self.plot_manager = PlotManagerFrame(self, self.controller)

        self.top_bar = SelectedExpsTopBar(self, self.controller)
        self.top_bar.pack(side="top", fill="x")

        self.exp_list_frame = ExpListFrame(self, self.controller, "single", self.on_select)
        self.exp_list_frame.config(padx=10, pady=10)
        self.exp_list_frame.pack(side="left", fill="y", pady=10)

        self.plot_manager.pack(side="right", fill="both", expand=True)

    def update_exp_list(self):
        mapped_exps = list(map(lambda exp: exp.title, self.controller.qram_exps))
        self.exp_list_frame.update_exp_list(mapped_exps)

    def on_select(self, event):
        idx = event.widget.curselection()[0]
        exp = self.controller.qram_exps[idx]
        self.open_popup = ChoosePlotPopup(self.controller, exp, self.plot_manager.selected_frame)
        self.open_popup.protocol("WM_DELETE_WINDOW", self.close_popup)

    def close_popup(self):
        self.open_popup.destroy()
        self.open_popup = None


class SelectedExpsTopBar(tk.Frame):
    def __init__(self, master, controller):
        super().__init__(master, padx=10, pady=10, bg="gray")
        self.master = master
        self.controller = controller

        tk.Button(self, text="<- Back", command=self.controller.back_to_search).pack(side="left")
        self.plus_button = tk.Button(self, text="+", command=self.plus_button_handler, padx=10)
        self.plus_button.pack(side="right")
        self.minus_button = tk.Button(self, text="-", command=self.minus_button_handler, padx=10, state="disabled")
        self.minus_button.pack(side="right")

    def plus_button_handler(self):
        disable = self.master.plot_manager.add_plot_frame()
        state = "disabled" if disable else "active"
        self.plus_button.config(state=state)
        self.minus_button.config(state="active")
    
    def minus_button_handler(self):
        disable = self.master.plot_manager.remove_plot_frame()
        state = "disabled" if disable else "active"
        self.minus_button.config(state=state)
        self.plus_button.config(state="active")



        