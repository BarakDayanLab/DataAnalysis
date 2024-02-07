import tkinter as tk
from data_visualization import QRAMDataVisualization

class ChoosePlotPopup(tk.Toplevel):
    def __init__(self, controller, exp, selected_frame):
        super().__init__(controller, takefocus=True)
        self.controller = controller
        self.exp = exp
        self.selected_frame = selected_frame

        self.geometry("800x400")
        self.title(f"Experiment time: {exp.title}")
        self.transient(controller)
        self.grab_set()

        func_dict = [func_name for func_name in dir(QRAMDataVisualization) if func_name.startswith("plot_")]

        self.selected_frame.current_exp = exp
        keys = self.selected_frame.plots.get(self.exp.title, {}).keys()
        self.buttons = {}
        for idx, func_name in enumerate(func_dict):
            row = idx // 4
            col = idx % 4
            txt = "del " if func_name in keys else ""
            self.buttons[func_name] = tk.Button(self, text=txt+func_name, command=lambda: self.button_pressed(func_name))
            self.buttons[func_name].grid(row=row, column=col)

    def button_pressed(self, func_name):
        has_plot = self.selected_frame.plots.get(self.exp.title, {}).get(func_name, None)
        func = getattr(self.selected_frame, func_name)
        if has_plot is None:
            self.buttons[func_name].config(text="del "+func_name)
            func()
        else:
            self.buttons[func_name].config(text=func_name)
            func(remove=True)
        


    # def plot_n_folded(self):
    #     self.plot_n_folded_button.config(text="Remove N Folded Seq", command=self.remove_n_folded)
    #     self.n_folded_plot = self.selected_frame.plot_n_folded(self.exp)

    # def remove_n_folded(self):
    #     self.n_folded_plot.pop().remove()
    #     self.plot_frames[-1].update()
    #     self.plot_n_folded_button.config(text="N Folded Seq", command=self.plot_n_folded)