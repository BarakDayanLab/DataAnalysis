import numpy as np
import tkinter as tk
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk


class QRAMDataVisualization:
    def plot_n_folded(self, remove=False):
        plot_id = "plot_n_folded"
        if remove:
            print("remove n folded")
            self.remove_plot(plot_id)
            return
        data = self.current_exp.files["North_timetags_folded_to_seq"] \
            + self.current_exp.files["Bright_timetags_folded_to_seq"] \
                + self.current_exp.files["Dark_timetags_folded_to_seq"] \
                    + self.current_exp.files["South_timetags_folded_to_seq"] \
                    +self.current_exp.files["FS_timetags_folded_to_seq"] 
        self.plot(plot_id, "North", "Counts", "Time (ns)", data, label=self.current_exp.title)


class PlotManagerFrame(tk.Frame):
    def __init__(self, master, controller):
        super().__init__(master)
        self.master = master
        self.controller = controller

        self.grid_rows = 2
        self.grid_columns = 3

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.plot_frames = [PlotFrame(self, self.controller, 0)]
        self.plot_frames[0].grid(row=0, column=0, sticky="nsew")
        self.selected_frame = self.plot_frames[0]


    def add_plot_frame(self):
        idx = len(self.plot_frames)
        row = idx // self.grid_columns
        col = idx % self.grid_columns

        self.plot_frames.append(PlotFrame(self, self.controller, idx))

        if col == 0:
            self.grid_rowconfigure(row, weight=1)
        if row == 0:
            self.grid_columnconfigure(col, weight=1)

        self.plot_frames[-1].grid(row=row, column=col, sticky="nsew")
        disable_button = True if idx == 5 else False
        return disable_button

    def remove_plot_frame(self):
        idx = len(self.plot_frames) - 1
        
        row = idx // self.grid_columns
        col = idx % self.grid_columns

        if row > (idx-1)//self.grid_columns and col == 0:
            self.grid_rowconfigure(row, weight=0)
        if col > (idx-1)%self.grid_columns and row == 0:
            self.grid_columnconfigure(col, weight=0)

        if self.selected_frame is self.plot_frames[-1]:
            self.selected_frame = self.plot_frames[-2]
        self.plot_frames.pop().destroy()

        disable_button = True if idx == 1 else False
        return disable_button
    
    def plot_frame_click_handler(self, plot_idx):
        self.selected_frame.config(highlightcolor="white")
        self.selected_frame = self.plot_frames[plot_idx]
        self.selected_frame.config(highlightcolor="black")



class PlotFrame(QRAMDataVisualization, tk.Frame):
    def __init__(self, master, controller, plot_id=None):
        super().__init__(master, highlightthickness=1, highlightcolor="white", padx=10, pady=10)
        self.master = master
        self.controller = controller
        self.id = plot_id
        self.plots = {}
        self.current_exp = None

        self.figure = Figure(dpi=100)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)

        self.widget = self.canvas.get_tk_widget()
        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.widget.pack(fill="both", expand=True)
        self.canvas.mpl_connect("button_release_event", lambda e: self.master.plot_frame_click_handler(self.id))
        
        self.ax = self.figure.add_subplot()
        self.ax.set_title(f"plot {self.id}")
        self.ax.set_ylabel("y")
        self.ax.set_xlabel("x")
        self.ax.legend() 

    def plot(self, plot_id, title=None, y_label="y", x_label="x", *args, **kwargs):
        self.ax.set_title(title or f"plot {self.id}")
        self.ax.set_ylabel(y_label)
        self.ax.set_xlabel(x_label)
        new_plot, = self.ax.plot(*args, **kwargs)

        exp_dict = self.plots.get(self.current_exp.title)

        if exp_dict is None:
            self.plots[self.current_exp.title] = {plot_id: new_plot}

        elif (old_plot:=exp_dict.get(plot_id)) is not None:
            old_plot.remove()
            self.plots[self.current_exp.title][plot_id] = new_plot
        else:
            self.plots[self.current_exp.title][plot_id] = new_plot

        self.ax.legend()
        self.canvas.draw()
    
    def remove_plot(self, plot_id):
        self.plots[self.current_exp.title].pop(plot_id).remove()
        self.ax.legend()
        self.canvas.draw()