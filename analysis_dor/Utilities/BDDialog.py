import tkinter as tk
from tkinter import filedialog


class BDDialog:

    def __init__(self, dialog_config=None):
        self.dialog_config = dialog_config
        self.res_button = None
        self.res_text = ''

    def set_dialog_config(self, dialog_config):
        self.dialog_config = dialog_config

    def handle_button(self, btn):

        # Remember the button id that was pressed
        self.res_button = btn

        if 'explorer' in btn and btn['explorer']:
            # Prevent an empty tkinter window from appearing
            #tk.Tk().withdraw()
            self.res_button['folder_name'] = filedialog.askdirectory()

        # Get the text from the text field and trim newlines and extra spaces at the start/end
        self.res_text = self.entry.get("1.0", tk.END)
        self.res_text = ' '.join(self.res_text.split())

        # Destroy the dialog window
        self.root.destroy()
        pass

    def prompt(self, dialog_config):

        title = dialog_config['title'] if 'title' in dialog_config else None
        message = dialog_config['message'] if 'message' in dialog_config else None
        geometry = dialog_config['geometry'] if 'geometry' in dialog_config else "500x270"

        # Create a root window, set title and geometry
        self.root = tk.Tk()
        if title is not None:
            self.root.title(title)
        self.root.geometry(geometry)
        self.root.wm_resizable(False, False)

        # Create a dialog window
        self.dialog = tk.Frame(self.root)

        self.dialog.pack(fill="both", expand=True)

        # Create a label with a message
        label = tk.Label(self.dialog, text=message, font=("Bold", 12))  # "Bold", 12)
        label.pack(pady=10)

        # Create the text field
        self.entry = tk.Text(self.dialog, height=10)
        self.entry.pack(side='top', padx=10, pady=10)

        # Create the buttons
        for btn in dialog_config['buttons']:
            padx = btn['pad'] if 'pad' in btn else 5
            button = tk.Button(self.dialog, text=btn['button_name'], font=("Bold", 12), command=lambda b=btn: self.handle_button(b))

            button.pack(side='left', padx=padx, pady=5)

        # Start the main loop
        self.root.mainloop()

        return self.res_text, self.res_button


if __name__ == "__main__":

    dialog_config = {
        "title": "This is a test title!",
        "geometry": "500x270",
        "message": "Enter the comment to save with experiment",
        "buttons": [
            {
                "description": "Saves results the are required for analysis",
                "button_id": 1,
                "button_name": "Analysis",
                "folder_name": "for_analysis"
            },
            {
                "description": "Saves results the are required for optimization",
                "button_id": 2,
                "button_name": "Optimization",
                "folder_name": "for_optimization"
            },
            {
                "description": "Ignore results - dont save",
                "button_id": 3,
                "button_name": "Ignore",
                "folder_name": "n/a"
            },
            {
                "description": "Allows generic save location",
                "button_id": 9,
                "button_name": "Choose Folder...",
                "pad": 40,
                "explorer": True
            }
        ]
    }

    bdd = BDDialog(dialog_config)
    res_text, res_button = bdd.prompt(dialog_config)

    pass
