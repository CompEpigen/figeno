import os
import customtkinter as ctk
import tkinter
import tkinter.messagebox
from PIL import Image
import importlib_resources as resources
from figeno.GUI.gui_colors import ColorButton, color2hovercolor
from figeno.GUI.gui_regions import RegionsFrame, HighlightsFrame
from figeno.GUI.gui_tracks import TracksFrame
from figeno.GUI.gui_utils import EntryPath
from figeno.tracks_plot import tracks_plot
import figeno.GUI.data
import json
font_title=("inter",16,"bold")
font=("inter",16)
font_small=("inter",10)

#TODO conda install -c conda-forge tk=*=xft_*

colorE=("gray60", "gray35")
colorD=("gray70", "gray25")
colorC=("gray85", "gray16")
colorB=("gray90", "gray13")
colorA=("gray95", "gray10")

colorE="gray35"
colorD="gray25"
colorC="gray16"
colorB="gray13"
colorA="gray10"

"""
colorA="#fcf5e0"
colorB="#fbf1d5"
colorC="#f9ebc4"
colorD="#f6e1a7"
colorE="#f1d892"
colorF="#ebcf81"
colorG="#e5c774"

BcolorA="#e2af1e"
BcolorB="#6c400c"
BcolorC="#6c400c"
"""


class GeneralFrame(ctk.CTkFrame):
    def __init__(self, master,params={}):
        super().__init__(master,fg_color=colorA,corner_radius=3,border_width=1)

        ctk.CTkLabel(self,10,20,text="General",font=font_title).grid(row=0,column=0,rowspan=1,padx=10,pady=(5,0),sticky="w")

        # Row 1: reference genome
        self.row1=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row1,10,20,text="Figure type: ",font=font).grid(row=0,column=0,padx=0,pady=5)
        self.type_menu=ctk.CTkOptionMenu(self.row1,values=["1 row", "2 rows", "circle"],font=font,dropdown_font=font,height=14,width=110)
        self.type_menu.grid(row=0,column=1,padx=(1,10),pady=5)
        ctk.CTkLabel(self.row1,10,20,text="Reference: ",font=font).grid(row=0,column=2,padx=0,pady=5)
        self.reference="hg19"
        self.reference_menu=ctk.CTkOptionMenu(self.row1,values=["hg19", "hg38", "custom"],
                                         command=self.optionmenu_reference,font=font,dropdown_font=font,height=14,width=110)
        self.reference_menu.grid(row=0,column=3,padx=(1,5),pady=5)
        self.chrarms_entry = EntryPath(self.row1,label="Chr arms file:",entry_width=120)
        self.genes_entry = EntryPath(self.row1,label="Genes file:",entry_width=120)
        if "chrarms_file" in params: self.chrarms_entry.set(params["chrarms_file"])
        if "genes_file" in params: self.genes_entry.set(params["genes_file"])
        if "reference" in params: 
            self.reference_menu.set(params["reference"])
            self.optionmenu_reference(params["reference"])
        
        # Row2: type
        #self.row2=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        #ctk.CTkLabel(self.row2,10,20,text="Figure type: ",font=font).grid(row=0,column=0,padx=0,pady=5)
        #self.type_menu=ctk.CTkOptionMenu(self.row2,values=["Linear", "Circular"],font=font,dropdown_font=font,height=14,width=110)
        #self.type_menu.grid(row=0,column=1,padx=(1,5),pady=5)

        self.row1.grid(row=1,column=0,columnspan=100,padx=10,pady=2,sticky="we")
        #self.row2.grid(row=2,column=0,columnspan=100,padx=10,pady=2,sticky="we")
        self.grid_columnconfigure(0,weight=0)
        self.grid_columnconfigure(1,weight=0)
        self.grid_columnconfigure(2,weight=0)

    def optionmenu_reference(self,choice):
        if choice !=self.reference:
            self.reference = choice
            if choice=="custom":
                self.chrarms_entry.grid(row=0,column=5,padx=(5,5),pady=5)
                self.genes_entry.grid(row=0,column=6,padx=(5,5),pady=5)
            else:
                self.chrarms_entry.grid_forget()
                self.genes_entry.grid_forget()  

    def reload(self,params={}):
        if "reference" in params:
            self.reference_menu.set(params["reference"])
            self.optionmenu_reference(params["reference"])
        else:
            self.reference_menu.set("hg19")
            self.optionmenu_reference("hg19")
        if "chrarms_file" in params: self.chrarms_entry.set(params["chrarms_file"])
        else: self.chrarms_entry.set("")
        if "genes_file" in params: self.genes_entry.set(params["genes_file"])
        else: self.genes_entry.set("")
        if "figure_type" in params: self.type_menu.set(params["figure_type"])
        else: self.type_menu.set("1 row")
        

        
    def get_params(self):
        params={"reference":self.reference_menu.get()}
        if self.reference_menu.get()=="custom":
            params["chrarms_file"] = self.chrarms_entry.get()
            params["genes_file"] = self.genes_entry.get()
        params["figure_type"] = self.type_menu.get()

        return params
    



class OutputFrame(ctk.CTkFrame):
    def __init__(self, master,params={}):
        super().__init__(master,fg_color=colorA,corner_radius=3,border_width=1)

        ctk.CTkLabel(self,10,20,text="Output",font=font_title).grid(row=0,column=0,rowspan=1,padx=10,pady=(5,0),sticky="w")
        self.row1=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row1,10,20,text="File: ",font=font).grid(row=0,column=0,padx=0,pady=5)
        self.path_stringvariable = tkinter.StringVar(self, "")
        if "file" in params: self.path_stringvariable.set(params["file"])
        self.entry_path=ctk.CTkEntry(self.row1,width=300,textvariable=self.path_stringvariable,font=font_small,height=12)
        self.entry_path.grid(row=0,column=1,padx=5,pady=5)
        with resources.as_file(resources.files(figeno.GUI.data) / "open.png") as infile:
            image_open = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,13))
        self.button_path = ctk.CTkButton(self.row1,text="",command=self.select_path,font=font,height=12,width=20,image=image_open)
        #self.button_path = ctk.CTkButton(self.row1,text="Browse",command=self.select_path,font=font,height=12,width=50)
        self.button_path.grid(row=0,column=2,padx=5,pady=5)
        ctk.CTkLabel(self.row1,10,20,text="dpi: ",font=font).grid(row=0,column=3,padx=(10,0),pady=5)
        self.dpi_stringvar = tkinter.StringVar(self, "400")
        if "dpi" in params: self.dpi_stringvar.set(str(params["dpi"]))
        self.dpi_entry=ctk.CTkEntry(self.row1,60,20,font=font,textvariable=self.dpi_stringvar)
        self.dpi_entry.grid(row=0,column=4,padx=(0,10))
        # Width
        ctk.CTkLabel(self.row1,10,20,text="Width (mm): ",font=font).grid(row=0,column=5,padx=(10,0),pady=5)
        self.width_stringvar = tkinter.StringVar(self, "183")
        if "width" in params: self.width_stringvar.set(str(params["width"]))
        self.width_entry=ctk.CTkEntry(self.row1,60,20,font=font,textvariable=self.width_stringvar)
        self.width_entry.grid(row=0,column=6,padx=(0,10))

        self.row1.grid(row=1,column=0,columnspan=100,padx=10,pady=5)

    def select_path(self):
        result = tkinter.filedialog.asksaveasfilename(initialfile="figure.svg")
        if len(result)>0:
            self.path_stringvariable.set(result)
            self.entry_path.xview(tkinter.END)

    def reload(self,params={}):
        if "file" in params: self.path_stringvariable.set(params["file"])
        else: self.path_stringvariable.set("")
        if "dpi" in params: self.dpi_stringvar.set(str(params["dpi"]))
        else: self.dpi_stringvar.set("400")
        if "width" in params: self.width_stringvar.set(str(params["width"]))
        else: self.width_stringvar.set("183")
    def get_params(self):
        params={"file":self.entry_path.get()}
        params["dpi"] = int(float(self.dpi_entry.get()))
        params["width"] = float(self.width_entry.get())

        return params
    


class App(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("figeno")
        self.geometry("1080x800")
        self.minsize(100,100)
        #self.maxsize(800,800)
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self.n_region_frames=1
        self.region_frames=[]
        self.config_file=None

        self.row1 = ctk.CTkFrame(self,fg_color="#111111",corner_radius=0)
        self.row1.grid(row=0,column=0,columnspan=100,sticky="nwe")
        with resources.as_file(resources.files(figeno.GUI.data) / "save.png") as infile:
            image_save = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(16,18))
        self.savebutton = ctk.CTkButton(self.row1,text="Save config",command=self.save_config,image=image_save,height=22,font=font)
        self.savebutton.grid(row=0,column=0,sticky="nw",padx=5,pady=5)
        with resources.as_file(resources.files(figeno.GUI.data) / "open.png") as infile:
            image_open = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(17,15))
        self.loadbutton = ctk.CTkButton(self.row1,text="Load config",command=self.load_config,image=image_open,height=22,font=font)
        self.loadbutton.grid(row=0,column=1,sticky="nw",padx=5,pady=5)
        with resources.as_file(resources.files(figeno.GUI.data) / "plus.png") as infile:
            image_plus = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(17,17))
        self.newbutton = ctk.CTkButton(self.row1,text="New",command=self.new_config,image=image_plus,height=22,font=font)
        self.newbutton.grid(row=0,column=2,sticky="nw",padx=5,pady=5)
        with resources.as_file(resources.files(figeno.GUI.data) / "run.png") as infile:
            image_run = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(17,17))
        self.runbutton = ctk.CTkButton(self.row1,text="Generate figure",image=image_run,fg_color="#016115",hover_color="#00420e",command=self.run,height=22,font=font)
        self.runbutton.grid(row=0,column=3,sticky="nw",padx=5,pady=5)



        self.main_frame = ctk.CTkScrollableFrame(self,fg_color=colorB,corner_radius=0)
        self.main_frame.grid(row=1,column=0,sticky="nwes",columnspan=100,pady=(0,5)) # ,sticky="news"

        self.general_frame = GeneralFrame(self.main_frame)
        self.general_frame.grid(row=0,column=0,padx=5,pady=4,sticky="we",columnspan=5)

        self.output_frame = OutputFrame(self.main_frame)
        self.output_frame.grid(row=1,column=0,padx=5,pady=4,sticky="we",columnspan=5)

        self.regions_frame = RegionsFrame(self.main_frame)
        self.regions_frame.grid(row=2,column=0,padx=5,pady=4,sticky="we",columnspan=5)

        self.highlights_frame = HighlightsFrame(self.main_frame)
        self.highlights_frame.grid(row=3,column=0,padx=5,pady=4,sticky="we",columnspan=5)

        self.tracks_frame = TracksFrame(self.main_frame)
        self.tracks_frame.grid(row=4,column=0,padx=5,pady=4,sticky="we",columnspan=5)
        self.main_frame.grid_columnconfigure(0,weight=1)

        

        self.rowconfigure(0,weight=0)
        self.rowconfigure(1,weight=1)
    
    def get_config(self):
        config={}
        config["general"] = self.general_frame.get_params()
        config["output"] = self.output_frame.get_params()
        config["regions"] = self.regions_frame.get_regions()
        config["highlights"] = self.highlights_frame.get_regions()
        config["tracks"] = self.tracks_frame.get_params_list()
        return config

    def save_config(self):
        config = self.get_config()
        if self.config_file is not None:
            filename = tkinter.filedialog.asksaveasfilename(initialfile=os.path.basename(self.config_file),initialdir=os.path.dirname(self.config_file))
        else:
            filename = tkinter.filedialog.asksaveasfilename(initialfile="config.json")
        if filename!="" and filename!="()":
            self.config_file = filename
            with open(filename,"w") as fp:
                json.dump(config,fp,indent= "\t")

    def load_config(self):
        filename = tkinter.filedialog.askopenfilename(filetypes=[("json","*.json")])
        if len(filename)>0:
            self.config_file=filename
            with open(filename,"r") as fp:
                config = json.load(fp)
                self.general_frame.reload(config["general"])
                self.output_frame.reload(config["output"])
                self.regions_frame.reload(config["regions"])
                if "highlights" in config: self.highlights_frame.reload(config["highlights"])
                else:  self.highlights_frame.reload()
                self.tracks_frame.reload(config["tracks"])
    def new_config(self):
        self.general_frame.reload()
        self.output_frame.reload()
        self.regions_frame.reload()
        self.highlights_frame.reload()
        self.tracks_frame.reload()
       
    def run(self):
        tp = tracks_plot(self.get_config())
        tp.draw()
        try:
            #tp = tracks_plot(self.get_config())
            tkinter.messagebox.showinfo("Success", "Figure was successfully generated.") 
        except Exception as e:
            tkinter.messagebox.showerror("Error", e) 
            print("Error during the figure generation.")
            print(e)



def main():
    ctk.set_appearance_mode("dark")
    app = App()
    app.mainloop()


if __name__=="__main__":
    main()



