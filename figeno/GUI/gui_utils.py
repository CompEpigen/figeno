import customtkinter as ctk
import tkinter
from PIL import Image, ImageTk
import importlib_resources as resources
import figeno.GUI.data
font=("inter",16)
font_small=("inter",10)


class EntryPath(ctk.CTkFrame):
    def __init__(self, master,label="",default_value="",entry_width=50,height=12):
        super().__init__(master,fg_color="transparent")
        if len(label)>0:
            ctk.CTkLabel(self,text=label,font=font).grid(row=0,column=0,padx=(0,2))
        self.stringvar = tkinter.StringVar(self, default_value)
        self.entry=ctk.CTkEntry(self,width=entry_width,textvariable=self.stringvar,font=font_small,height=height)
        self.entry.grid(row=0,column=1)
        with resources.as_file(resources.files(figeno.GUI.data) / "open.png") as infile:
            image_open = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,13))
        self.button_gtf = ctk.CTkButton(self,text="",image=image_open,command=self.select_path,font=font,height=12,width=20)
        self.button_gtf.grid(row=0,column=2,padx=2)

    def select_path(self):
        result = tkinter.filedialog.askopenfilename()
        if len(result)>0:
            self.stringvar.set(result)
            self.entry.xview(tkinter.END)

    def get(self):
        return self.stringvar.get()
    def set(self,val):
        self.stringvar.set(val)
        self.entry.xview(tkinter.END)