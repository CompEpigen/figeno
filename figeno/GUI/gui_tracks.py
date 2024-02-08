import customtkinter as ctk
import tkinter
from PIL import Image
import importlib_resources as resources
from figeno.GUI.gui_colors import ColorSelector, ColorButton, color2hovercolor
from figeno.GUI.gui_utils import EntryPath
import figeno.GUI.data
font_title=("inter",16,"bold")
font=("inter",16)
font_small=("inter",10)


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


class TracksFrame(ctk.CTkFrame):
    def __init__(self, master,params_list=[]):
        super().__init__(master,width=10,fg_color=colorA,corner_radius=3,border_width=1)
        self.track_frames=[]
        self.unused_track_frames=[]

        self.label = ctk.CTkLabel(self,text="Tracks",font=font_title)
        self.label.grid(row=0,column=0,padx=10,pady=5,sticky="w")
        with resources.as_file(resources.files(figeno.GUI.data) / "plus.png") as infile:
            image_plus = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(17,17))
        self.button =ctk.CTkButton(self, text="Add track",font=font, command=self.add_track_frame,image=image_plus,height=20)
        self.button.grid(row=0, column=1, padx=10, pady=5, sticky="w")
        self.columnconfigure((0,1,2), weight=0)
        self.columnconfigure(2, weight=1)
        self.button2 =ctk.CTkButton(self, text="Load files",font=font, command=self.load_files,height=20)
        self.button2.grid(row=0, column=2, padx=10, pady=5, sticky="w")
        for params in params_list[::-1]:
            self.add_track_frame(params)

    def add_track_frame(self,params={},update_grid=True,append=False):
        if len(self.unused_track_frames)>0:
            track_frame = self.unused_track_frames.pop()
            track_frame.reload(params)
        else:
            track_frame = TrackFrame(self,params)
        if append:
            self.track_frames.append(track_frame)
        else:
            self.track_frames = [track_frame] + self.track_frames
        if update_grid:
            self.update_trackframes_grid()

    def remove_track_frame(self,trackframe):
        trackframe.grid_forget()
        self.unused_track_frames.append(trackframe)
        self.track_frames.remove(trackframe)

    def move_frame_up(self,frame):
        for i in range(1,len(self.track_frames)):
            if self.track_frames[i]==frame:
                self.track_frames[i-1],self.track_frames[i] = self.track_frames[i],self.track_frames[i-1]
                break
        self.update_trackframes_grid()

    def move_frame_down(self,frame):
        for i in range(len(self.track_frames)-1):
            if self.track_frames[i]==frame:
                self.track_frames[i+1],self.track_frames[i] = self.track_frames[i],self.track_frames[i+1]
                break
        self.update_trackframes_grid()
    def move_frame_upup(self,frame):
        self.track_frames.remove(frame)
        self.track_frames = [frame] + self.track_frames
        self.update_trackframes_grid()

    def move_frame_downdown(self,frame):
        self.track_frames.remove(frame)
        self.track_frames =  self.track_frames + [frame]
        self.update_trackframes_grid()

    def update_trackframes_grid(self):
        for i in range(len(self.track_frames)):
            pady=(4,4)
            if i==0: pady =  (5,pady[1])
            if i==len(self.track_frames)-1: pady =  (pady[0],5)
            self.track_frames[i].grid(row=i+1, column=0, padx=10, pady=pady, sticky="new",columnspan=3)

    def load_files(self):
        result = tkinter.filedialog.askopenfilenames()
        if len(result)>0:
            for file in result:
                if file.endswith(".bam"): self.add_track_frame({"type":"alignments","bam":file})
                elif file.endswith(".bw") or file.endswith(".bigwig") or file.endswith(".bigWig"):  self.add_track_frame({"type":"bigwig","bigwig":file})
                elif file.endswith("cool"): self.add_track_frame({"type":"HiC","file":file})
                elif file.endswith(".gtf") or file.endswith(".gtf.gz"): self.add_track_frame({"type":"genes","gtf":file})
                elif file.endswith(".bed"): pass # TODO
                
    def reload(self,params_list=[]):
        n_previous_tracks = len(self.track_frames)
        n_new_tracks = len(params_list)
        if n_previous_tracks>n_new_tracks:
            for i in range(n_previous_tracks-n_new_tracks):
                track_frame = self.track_frames.pop()
                track_frame.grid_forget()
                self.unused_track_frames.append(track_frame)
        for i in range(min(n_new_tracks,n_previous_tracks)):
            self.track_frames[i].reload(params_list[i])
        if n_new_tracks>n_previous_tracks:
            for i in range(n_previous_tracks,n_new_tracks):
                self.add_track_frame(params_list[i],update_grid=False,append=True)
        self.update_trackframes_grid()
            
    def print_tracks(self):
        for x in self.tracks_frames:
            print(x.get_track())
    def get_params_list(self):
        l=[]
        for x in self.track_frames: l.append(x.get_params())
        return l

class TrackFrame(ctk.CTkFrame):
    def __init__(self, master,params={}):
        super().__init__(master,fg_color=colorB,border_width=2)

        self.options_frame=None
        self.current_type = " "
        with resources.as_file(resources.files(figeno.GUI.data) / "up.png") as infile:
            image_up = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,15))
        with resources.as_file(resources.files(figeno.GUI.data) / "down.png") as infile:
            image_down = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,15))
        self.up_button=ctk.CTkButton(self,width=5,height=5,text="",image=image_up,command=self.up,font=font,fg_color=colorA,hover_color=colorC)
        self.down_button=ctk.CTkButton(self,5,5,text="",image=image_down,command=self.down,font=font,fg_color=colorA,hover_color=colorC)

        with resources.as_file(resources.files(figeno.GUI.data) / "upup.png") as infile:
            image_upup =ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,15))
        with resources.as_file(resources.files(figeno.GUI.data) / "downdown.png") as infile:
            image_downdown = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,15))
        self.upup_button=ctk.CTkButton(self,width=5,height=5,text="",image=image_upup,command=self.upup,font=font,fg_color=colorA,hover_color=colorC)
        self.downdown_button=ctk.CTkButton(self,5,5,text="",image=image_downdown,command=self.downdown,font=font,fg_color=colorA,hover_color=colorC)

        typelabel = ctk.CTkLabel(self,text="Track type:",font=font)
        #typevar = ctk.StringVar(value=self.current_type)
        self.typemenu = ctk.CTkOptionMenu(self,values=["", "chr axis", "genes" ,"bed" ,"bigwig", "coverage", "alignments", "Met freq", "HiC","SV", "copynumber", "ase"],
                                         command=self.optionmenu_type,font=font,dropdown_font=font_small,height=14,width=110)
        
        self.height_label = ctk.CTkLabel(self,text="Height (mm):",font=font)
        self.height_stringvariable = tkinter.StringVar(self, "10")
        self.height_entry=ctk.CTkEntry(self,width=60,textvariable=self.height_stringvariable,font=font,height=12)

        self.margin_label = ctk.CTkLabel(self,text="Margin above:",font=font)
        self.margin_stringvariable = tkinter.StringVar(self, "1.5")
        if "margin_above" in params: self.margin_stringvariable.set(str(params["margin_above"]))
        self.margin_entry=ctk.CTkEntry(self,width=40,textvariable=self.margin_stringvariable,font=font,height=12)

        self.boundingbox_label = ctk.CTkLabel(self,text="Box:",font=font)
        self.boundingbox_checkbox=ctk.CTkCheckBox(self,width=20,height=20,checkbox_width=20,checkbox_height=20,text="")

        self.label_fontscale=ctk.CTkLabel(self,10,10,text="Fontscale:",font=font)
        self.fontscale_stringvariable = tkinter.StringVar(self, "1.0")
        if "fontscale" in params: self.fontscale_stringvariable.set(str(params["fontscale"]))
        self.fontscale_entry = ctk.CTkEntry(self, width=60,height=20,font=font,textvariable=self.fontscale_stringvariable)

        with resources.as_file(resources.files(figeno.GUI.data) / "delete.png") as infile:
            image_delete =ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS))
        self.delete_button=ctk.CTkButton(self,10,10,text="",image=image_delete,command=self.removeframe,font=font,fg_color="#e74c3c",hover_color=color2hovercolor("#e74c3c")) # colorD #e74c3c #ff291e

        self.up_button.grid(row=0,column=0,sticky="n",padx=2,pady=2)
        self.down_button.grid(row=1,column=0,sticky="s",padx=2,pady=3)
        self.upup_button.grid(row=0,column=1,sticky="n",padx=2,pady=2)
        self.downdown_button.grid(row=1,column=1,sticky="s",padx=2,pady=3)
        typelabel.grid(row=0,column=2,padx=(10,2),pady=4)
        self.typemenu.grid(row=0,column=3,pady=4)

        self.delete_button.grid(row=0,column=20,rowspan=2,padx=10)

        self.columnconfigure((0,1,2,3,5,6,7,8,9,10,11,20), weight=0)
        self.columnconfigure(12, weight=1)
        if "type" in params:
            self.optionmenu_type(params["type"],params)
            self.current_type=params["type"]
            self.typemenu.set(params["type"])

        
    def removeframe(self):
        self.master.remove_track_frame(self)
    def up(self):
        self.master.move_frame_up(self)
    def down(self):
        self.master.move_frame_down(self)
    def upup(self):
        self.master.move_frame_upup(self)
    def downdown(self):
        self.master.move_frame_downdown(self)

    def optionmenu_type(self,value,params={}):
        if self.current_type!=value:
            self.current_type=value
            if self.options_frame is not None:
                self.options_frame.destroy()

            if value=="genes":
                self.options_frame=GenesOptionsFrame(self,params)
            elif value=="bed":
                self.options_frame=BedOptionsFrame(self,params)
            elif value=="bigwig":
                self.options_frame=BigwigOptionsFrame(self,params)
            elif value=="coverage":
                self.options_frame=CoverageOptionsFrame(self,params)
            elif value=="chr axis":
                self.options_frame=ChrOptionsFrame(self,params)
            elif value=="alignments":
                self.options_frame=AlignmentsOptionsFrame(self,params)
            elif value=="Met freq":
                self.options_frame=MetfreqOptionsFrame(self,params)
            elif value=="HiC":
                self.options_frame=HiCOptionsFrame(self,params)
            elif value=="SV":
                self.options_frame=SVOptionsFrame(self,params)
            elif value=="copynumber":
                self.options_frame=CopynumberOptionsFrame(self,params)
            elif value=="ase":
                self.options_frame=AseOptionsFrame(self,params)
            else:
                self.options_frame = None
            
            # Set height
            if "height" in params: self.height_stringvariable.set(str(params["height"]))
            else:
                if value in ["chr axis","bigwig","coverage","bed"]: self.height_stringvariable.set("10")
                elif value in ["genes"]:  self.height_stringvariable.set("12")
                elif value in ["HiC"]:  self.height_stringvariable.set("50")
                elif value in ["alignments"]:  self.height_stringvariable.set("75")
                elif value in ["Met freq"]:  self.height_stringvariable.set("20")
                elif value in ["copynumber"]: self.height_stringvariable.set("30")
                elif value in ["SV"]: self.height_stringvariable.set("15")

            # Set bounding box
            if "bounding_box" in params:
                if params["bounding_box"]: self.boundingbox_checkbox.select()
                else: self.boundingbox_checkbox.deselect()
            else:
                if value in ["alignments","Met freq","SV","copynumber"]: self.boundingbox_checkbox.select()
                else: self.boundingbox_checkbox.deselect()

            # Add elements to grid
            if value!="":
                self.height_label.grid(row=0,column=4,pady=4,padx=(10,0))
                self.height_entry.grid(row=0,column=5,pady=4,padx=5)
                self.margin_label.grid(row=0,column=6,pady=4,padx=(10,0))
                self.margin_entry.grid(row=0,column=7,pady=4,padx=5)
                self.boundingbox_label.grid(row=0,column=8,pady=4,padx=5)
                self.boundingbox_checkbox.grid(row=0,column=9,pady=4,padx=0)
                self.label_fontscale.grid(row=0,column=10,pady=4,padx=0)
                self.fontscale_entry.grid(row=0,column=11,pady=4,padx=0)
                self.options_frame.grid(row=1,column=2,sticky="nesw",columnspan=12,padx=10,pady=4)
            else:
                self.height_label.grid_forget()
                self.height_entry.grid_forget()
                self.margin_label.grid_forget()
                self.margin_entry.grid_forget()
                self.boundingbox_label.grid_forget()
                self.boundingbox_checkbox.grid_forget()
                self.label_fontscale.grid_forget()
                self.fontscale_entry.grid_forget()
                if self.options_frame is not None:
                    self.options_frame.grid_forget()

    def reload(self,params={}):
        self.current_type = "reload"
        if "type" in params: 
            self.optionmenu_type(params["type"],params)
            self.typemenu.set(params["type"])
        else: 
            self.optionmenu_type("")
            self.typemenu.set("")


    def get_params(self):
        if self.typemenu.get()=="": return {}
        
        params = {}
        params["type"] = self.typemenu.get()
        params["height"] = float(self.height_stringvariable.get())
        params["margin_above"] = float(self.margin_stringvariable.get())
        params["bounding_box"] = self.boundingbox_checkbox.get() ==1
        params["fontscale"] = float(self.fontscale_entry.get())
        params_opt = self.options_frame.get_params()
        for k in params_opt:
            params[k] = params_opt[k]
        return params

    
class ChrOptionsFrame(ctk.CTkFrame):
    def __init__(self, master,params={}):
        super().__init__(master,fg_color="transparent",border_width=0)

        self.label_style=ctk.CTkLabel(self,10,10,text="Style: ",font=font)
        self.stylemenu = ctk.CTkOptionMenu(self,values=["Default", "Arrow", "Ideogram"],font=font,dropdown_font=font,height=12)
        if "style" in params: self.stylemenu.set(params["style"])
        
        self.label_style.grid(row=0,column=0,pady=4)
        self.stylemenu.grid(row=0,column=1,padx=10,pady=4)

        ctk.CTkLabel(self,10,10,text="Unit: ",font=font).grid(row=0,column=2,pady=4,padx=(5,0))
        self.unit_menu = ctk.CTkOptionMenu(self,values=["bp", "kb", "Mb"],font=font,dropdown_font=font,height=12,width=50)
        if "unit" in params: self.unit_menu.set(params["unit"])
        else: self.unit_menu.set("kb")
        self.unit_menu.grid(row=0,column=3,pady=4,padx=(1,5))
        ctk.CTkLabel(self,10,10,text="Ticks position: ",font=font).grid(row=0,column=4,pady=4,padx=(5,0))
        self.ticks_menu = ctk.CTkOptionMenu(self,values=["below", "above","none"],font=font,dropdown_font=font,height=12,width=80)
        if "ticks_pos" in params: self.ticks_menu.set(params["ticks_pos"])
        self.ticks_menu.grid(row=0,column=5,pady=4,padx=(1,5))


    def get_params(self):
        params =  {"style":self.stylemenu.get()}
        params["unit"] = self.unit_menu.get()
        params["ticks_pos"] = self.ticks_menu.get()
        return params
    
    

    
class GenesOptionsFrame(ctk.CTkFrame):
    def __init__(self, master,params={},**kwargs):
        super().__init__(master,fg_color="transparent",border_width=0)

        # Row1: style and exon color and gene _names
        self.row1=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        self.label_style=ctk.CTkLabel(self.row1,10,10,text="Style:",font=font)
        self.stylemenu = ctk.CTkOptionMenu(self.row1,values=["Default", "TSS_arrow"],font=font,dropdown_font=font,height=12)
        if "style" in params: self.stylemenu.set(params["style"])
        self.label_style.grid(row=0,column=0,pady=4,sticky="w")
        self.stylemenu.grid(row=0,column=1,padx=10,pady=4,sticky="w")
        ctk.CTkLabel(self.row1,10,10,text="Exon color:",font=font).grid(row=0,column=2,pady=4,padx=(5,0))

        self.buttoncolor = ColorButton(self.row1,25,25,color="#2178b0")
        if "exon_color" in params: self.buttoncolor.set(params["exon_color"])
        self.buttoncolor.grid(row=0,column=3,padx=5,pady=4)

        ctk.CTkLabel(self.row1,10,10,text="Genes: ",font=font).grid(row=0,column=4,pady=4,padx=(5,0))
        self.genes_stringvariable = tkinter.StringVar(self, "auto")
        if "gene_names" in params: self.genes_stringvariable.set(params["gene_names"])
        self.entry_genes=ctk.CTkEntry(self.row1,width=100,textvariable=self.genes_stringvariable,font=font,height=12)
        self.entry_genes.grid(row=0,column=5,pady=4,padx=(1,5))
        #self.colorRow = colorFrame(self,"Color")        
        
        self.row1.grid(row=0,column=0,pady=0,sticky="ew")
        #self.colorRow.grid(row=2,column=0,pady=0,sticky="ew")

    def get_params(self):
        params = {}
        params["style"] = self.stylemenu.get()
        params["exon_color"] = self.buttoncolor.get()
        params["gene_names"] = self.entry_genes.get()
        return params
    
class BedOptionsFrame(ctk.CTkFrame):
    def __init__(self, master,params={},**kwargs):
        super().__init__(master,fg_color="transparent",border_width=0)


        # Row1: file path
        self.row1=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        self.file_entry = EntryPath(self.row1,label="File: ",entry_width=300)
        if "bed" in params: self.file_entry.set(params["bed"])
        self.file_entry.grid(row=0,column=1,pady=5,padx=(0,5))

        # Row2: label
        self.row2=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row2,10,10,text="Label: ",font=font).grid(row=0,column=0,pady=5)
        self.label_stringvariable = tkinter.StringVar(self, "")
        if "label" in params: self.label_stringvariable.set(params["label"])
        self.label_entry=ctk.CTkEntry(self.row2,width=200,textvariable=self.label_stringvariable,font=font,height=12)
        self.label_entry.grid(row=0,column=1,padx=5,pady=5)
        ctk.CTkLabel(self.row2,10,10,text="Rotate label:",font=font).grid(row=0,column=2,padx=5,pady=5)
        self.label_switch = ctk.CTkSwitch(self.row2,10,10,text="",font=font)
        if "label_rotate" in params and params["label_rotate"]: self.label_switch.select()
        self.label_switch.grid(row=0,column=3,padx=5,pady=5)



        # Colors
        ctk.CTkLabel(self,10,10,text="Color",font=font).grid(row=0,column=1,padx=(20,5),pady=5)
        self.buttoncolor = ColorButton(self,40,40,color="#2980b9")
        if "color" in params: self.buttoncolor.set(params["color"])
        self.buttoncolor.grid(row=1,column=1,padx=(20,5))
        #self.colorRow = colorFrame(self,"Color")        
        
        self.row1.grid(row=0,column=0,pady=0,sticky="ew")
        self.row2.grid(row=1,column=0,pady=0,sticky="ew")

    def get_params(self):
        params =  {"bed":self.file_entry.get()}
        params["label"] = self.label_entry.get()
        params["label_rotate"] = self.label_switch.get() == 1
        params["color"] = self.buttoncolor.get()
        return params


class BigwigOptionsFrame(ctk.CTkFrame):
    def __init__(self, master,params={},**kwargs):
        super().__init__(master,fg_color="transparent",border_width=0)


        # Row1: file path and n_bins
        self.row1=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        self.file_entry = EntryPath(self.row1,label="File: ",entry_width=300)
        if "bigwig" in params: self.file_entry.set(params["bigwig"])
        self.file_entry.grid(row=0,column=1,pady=5,padx=(0,5))
        ctk.CTkLabel(self.row1,10,10,text="n_bins: ",font=font).grid(row=0,column=3,pady=5,padx=(15,0))
        self.nbins_stringvariable = tkinter.StringVar(self, "500")
        if "n_bins" in params: self.nbins_stringvariable.set(str(params["n_bins"]))
        self.n_bins_entry = ctk.CTkEntry(self.row1,width=80,textvariable=self.nbins_stringvariable,font=font,height=12)
        self.n_bins_entry.grid(row=0,column=4,padx=5,pady=5)

        # Row2: label
        self.row2=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row2,10,10,text="Label: ",font=font).grid(row=0,column=0,pady=5)
        self.label_stringvariable = tkinter.StringVar(self, "")
        if "label" in params: self.label_stringvariable.set(params["label"])
        self.label_entry=ctk.CTkEntry(self.row2,width=200,textvariable=self.label_stringvariable,font=font,height=12)
        self.label_entry.grid(row=0,column=1,padx=5,pady=5)
        ctk.CTkLabel(self.row2,10,10,text="Rotate label:",font=font).grid(row=0,column=2,padx=5,pady=5)
        self.label_switch = ctk.CTkSwitch(self.row2,10,10,text="",font=font)
        if "label_rotate" in params and params["label_rotate"]: self.label_switch.select()
        self.label_switch.grid(row=0,column=3,padx=5,pady=5)


        # Row3: scale
        self.row3=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row3,10,10,text="Scale: ",font=font).grid(row=0,column=0,pady=5)
        self.scale="auto"
        self.scalemenu = ctk.CTkOptionMenu(self.row3,values=["auto", "auto per region", "custom"],
                                         command=self.scalemenu_callback,font=font,dropdown_font=font,height=12)
        if "scale" in params:
            self.scale = params["scale"]
            self.scalemenu.set(params["scale"])
        self.scalemenu.grid(row=0,column=1,padx=5,pady=5)
        ctk.CTkLabel(self.row3,10,10,text="Max: ",font=font).grid(row=0,column=2,padx=5,pady=5)
        self.max_stringvariable = tkinter.StringVar(self, "auto")
        if "scale_max" in params and self.scale=="custom": self.max_stringvariable.set(str(params["scale_max"]))
        self.max_entry=ctk.CTkEntry(self.row3,width=50,textvariable=self.max_stringvariable,font=font,height=12,state="disabled")
        if "scale" in params and params["scale"]=="custom":self.max_entry.configure(state="normal")
        self.max_entry.grid(row=0,column=3,padx=(0,5),pady=5)
        ctk.CTkLabel(self.row3,10,10,text="Position: ",font=font).grid(row=0,column=4,padx=5,pady=5)
        self.scalepos_menu = ctk.CTkOptionMenu(self.row3,values=["left","corner","corner all","none"],font=font,dropdown_font=font,height=12,width=80)
        if "scale_pos" in params: self.scalepos_menu.set(params["scale_pos"])
        if self.scale=="auto per region": self.scalepos_menu.configure(values=["corner all"])
        self.scalepos_menu.grid(row=0,column=5,padx=(0,5),pady=5)

        # Colors
        ctk.CTkLabel(self,10,10,text="Color",font=font).grid(row=0,column=1,padx=(20,5),pady=5)
        self.buttoncolor = ColorButton(self,40,40,color="#2980b9")
        if "color" in params: self.buttoncolor.set(params["color"])
        self.buttoncolor.grid(row=1,column=1,padx=(20,5))
        #self.colorRow = colorFrame(self,"Color")        
        
        self.row1.grid(row=0,column=0,pady=0,sticky="ew")
        self.row2.grid(row=1,column=0,pady=0,sticky="ew")
        self.row3.grid(row=2,column=0,pady=0,sticky="ew")

    def scalemenu_callback(self,choice):
        if choice!=self.scale:
            if choice=="auto per region": 
                self.scalepos_menu.configure(values=["corner all","none"])
                if self.scalepos_menu.get()!="none":
                    self.scalepos_menu.set("corner all")
            elif self.scale=="auto per region": 
                self.scalepos_menu.configure(values=["left","corner","corner all","none"])
                self.scalepos_menu.set("left")
            
            if choice=="custom": 
                self.max_entry.configure(state="normal")
                self.max_stringvariable.set("1")
            else:
                self.max_entry.configure(state="disabled")
                self.max_stringvariable.set("auto")
            self.scale=choice
    def get_params(self):
        params =  {"bigwig":self.file_entry.get()}
        params["n_bins"] = int(self.n_bins_entry.get())
        params["label"] = self.label_entry.get()
        params["label_rotate"] = self.label_switch.get() == 1
        params["scale"] = self.scalemenu.get()
        if self.scalemenu.get()=="custom": 
            if "," in self.max_entry.get():
                params["scale_max"] = self.max_entry.get()
            else:
                params["scale_max"] = float(self.max_entry.get())
        params["scale_pos"] = self.scalepos_menu.get()
        params["color"] = self.buttoncolor.get()


        return params
    
class CoverageOptionsFrame(ctk.CTkFrame):
    def __init__(self, master,params={},**kwargs):
        super().__init__(master,fg_color="transparent",border_width=0)


        # Row1: file path and n_bins
        self.row1=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        self.file_entry = EntryPath(self.row1,label="File: ",entry_width=300)
        if "bam" in params: self.file_entry.set(params["bam"])
        self.file_entry.grid(row=0,column=0,pady=5,padx=(0,5))
        ctk.CTkLabel(self.row1,10,10,text="n_bins: ",font=font).grid(row=0,column=3,pady=5,padx=(15,0))
        self.nbins_stringvariable = tkinter.StringVar(self, "500")
        if "n_bins" in params: self.nbins_stringvariable.set(str(params["n_bins"]))
        self.n_bins_entry = ctk.CTkEntry(self.row1,width=80,textvariable=self.nbins_stringvariable,font=font,height=12)
        self.n_bins_entry.grid(row=0,column=4,padx=5,pady=5)

        # Row2: label
        self.row2=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row2,10,10,text="Label: ",font=font).grid(row=0,column=0,pady=5)
        self.label_stringvariable = tkinter.StringVar(self, "")
        if "label" in params: self.label_stringvariable.set(params["label"])
        self.label_entry=ctk.CTkEntry(self.row2,width=200,textvariable=self.label_stringvariable,font=font,height=12)
        self.label_entry.grid(row=0,column=1,padx=5,pady=5)
        ctk.CTkLabel(self.row2,10,10,text="Rotate label:",font=font).grid(row=0,column=2,padx=5,pady=5)
        self.label_switch = ctk.CTkSwitch(self.row2,10,10,text="",font=font)
        if "label_rotate" in params and params["label_rotate"]: self.label_switch.select()
        self.label_switch.grid(row=0,column=3,padx=5,pady=5)


        # Row3: scale
        self.row3=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row3,10,10,text="Scale: ",font=font).grid(row=0,column=0,pady=5)
        self.scale="auto"
        self.scalemenu = ctk.CTkOptionMenu(self.row3,values=["auto", "auto per region", "custom"],
                                         command=self.scalemenu_callback,font=font,dropdown_font=font,height=12)
        if "scale" in params:
            self.scale = params["scale"]
            self.scalemenu.set(params["scale"])
        self.scalemenu.grid(row=0,column=1,padx=5,pady=5)
        ctk.CTkLabel(self.row3,10,10,text="Max: ",font=font).grid(row=0,column=2,padx=5,pady=5)
        self.max_stringvariable = tkinter.StringVar(self, "auto")
        if "scale_max" in params and self.scale=="custom": self.max_stringvariable.set(str(params["scale_max"]))
        self.max_entry=ctk.CTkEntry(self.row3,width=50,textvariable=self.max_stringvariable,font=font,height=12,state="disabled")
        self.max_entry.grid(row=0,column=3,padx=(0,5),pady=5)
        ctk.CTkLabel(self.row3,10,10,text="Position: ",font=font).grid(row=0,column=4,padx=5,pady=5)
        self.scalepos_menu = ctk.CTkOptionMenu(self.row3,values=["left","corner","corner all","none"],font=font,dropdown_font=font,height=12)
        if "scale_pos" in params: self.scalepos_menu.set(params["scale_pos"])
        if self.scale=="auto per region": self.scalepos_menu.configure(values=["corner all","none"])
        self.scalepos_menu.grid(row=0,column=5,padx=(0,5),pady=5)

        # Colors
        ctk.CTkLabel(self,10,10,text="Color",font=font).grid(row=0,column=1,padx=5,pady=5)
        self.buttoncolor= ColorButton(self,40,40,color="#888888")
        self.buttoncolor.grid(row=1,column=1)
        #self.colorRow = colorFrame(self,"Color")        
        
        self.row1.grid(row=0,column=0,pady=0,sticky="ew")
        self.row2.grid(row=1,column=0,pady=0,sticky="ew")
        self.row3.grid(row=2,column=0,pady=0,sticky="ew")

    def scalemenu_callback(self,choice):
        if choice!=self.scale:
            if choice=="auto per region": 
                self.scalepos_menu.configure(values=["corner all","none"])
                self.scalepos_menu.set("corner all")
            elif self.scale=="auto per region": 
                self.scalepos_menu.configure(values=["left","corner","corner all","none"])
                self.scalepos_menu.set("left")
            
            if choice=="custom": 
                self.max_entry.configure(state="normal")
                self.max_stringvariable.set("1")
            else:
                self.max_entry.configure(state="disabled")
                self.max_stringvariable.set("auto")
            self.scale=choice
    def get_params(self):
        params =  {"bam":self.file_entry.get()}
        params["n_bins"] = int(self.n_bins_entry.get())
        params["label"] = self.label_entry.get()
        params["label_rotate"] = self.label_switch.get() == 1
        params["scale"] = self.scalemenu.get()
        if self.scalemenu.get=="custom":
            params["scale_max"] = float(self.max_entry.get())
        params["scale_pos"] = self.scalepos_menu.get()
        params["color"] = self.buttoncolor.get()


        return params



class AlignmentsOptionsFrame(ctk.CTkFrame):
    def __init__(self, master,params={},**kwargs):
        super().__init__(master,fg_color="transparent",border_width=0)


        # Row1: file path
        self.row1=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        self.file_entry = EntryPath(self.row1,label="File: ",entry_width=300)
        if "bam" in params: self.file_entry.set(params["bam"])
        self.file_entry.grid(row=0,column=0,padx=(0,5),pady=5)

        # Row2: label
        self.row2=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row2,10,10,text="Label: ",font=font).grid(row=0,column=0,pady=5)
        self.label_stringvariable = tkinter.StringVar(self, "")
        if "label" in params: self.label_stringvariable.set(params["label"])
        self.label_entry=ctk.CTkEntry(self.row2,width=200,textvariable=self.label_stringvariable,font=font,height=12)
        self.label_entry.grid(row=0,column=1,padx=5,pady=5)
        ctk.CTkLabel(self.row2,10,10,text="Rotate label:",font=font).grid(row=0,column=2,padx=5,pady=5)
        self.label_switch = ctk.CTkSwitch(self.row2,10,10,text="",font=font)
        if "label_rotate" in params and params["label_rotate"]: self.label_switch.select()
        self.label_switch.grid(row=0,column=3,padx=5,pady=5)

        # Row3: group by. Either nothing or haplotype (if haplotype, show unphased or not ?). Add other options ? SNPs ? Maybe allow using custom functions
        self.row3=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row3,10,10,text="Group by: ",font=font).grid(row=0,column=0,pady=5)
        self.group_by="none"
        self.groupby_menu = ctk.CTkOptionMenu(self.row3,values=["none", "haplotype"],font=font,dropdown_font=font,height=12,command=self.update_groupby)
        self.groupby_menu.grid(row=0,column=1,pady=5)

        self.exchangehaplotypes_label =  ctk.CTkLabel(self.row3,10,10,text="Exchange haplotypes: ",font=font)
        self.exchangehaplotypes_checkbox = ctk.CTkCheckBox(self.row3,width=20,height=20,checkbox_width=20,checkbox_height=20,text="")
        if "exchange_haplotypes" in params and params["exchange_haplotypes"]: self.exchangehaplotypes_checkbox.select()
        self.showunphased_label =  ctk.CTkLabel(self.row3,10,10,text="Show unphased: ",font=font)
        self.showunphased_checkbox = ctk.CTkCheckBox(self.row3,width=20,height=20,checkbox_width=20,checkbox_height=20,text="",command=self.update_colors_labels)
        if (not "show_unphased" in params) or params["show_unphased"]: self.showunphased_checkbox.select()
        self.colors_label =  ctk.CTkLabel(self.row3,10,10,text="Colors: ",font=font)
        self.colors_checkbox = ctk.CTkCheckBox(self.row3,width=20,height=20,checkbox_width=20,checkbox_height=20,text="",command=self.update_colors_labels)
        if (not "color_haplotypes" in params) or params["color_haplotypes"]: self.colors_checkbox.select()


        


        self.row_grouplabels=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row_grouplabels,10,10,text="Group labels: ",font=font).grid(row=0,column=0,pady=5)
        self.grouplabel1_stringvariable = tkinter.StringVar(self, "Haplotype 1")
        self.grouplabel2_stringvariable = tkinter.StringVar(self, "Haplotype 2")
        self.grouplabel3_stringvariable = tkinter.StringVar(self, "Unphased")
        if "group_labels" in params:
            self.grouplabel1_stringvariable.set(params["group_labels"][0])
            self.grouplabel2_stringvariable.set(params["group_labels"][1])
            if len(params["group_labels"])>2: self.grouplabel3_stringvariable.set(params["group_labels"][2])
        self.grouplabel1_entry=ctk.CTkEntry(self.row_grouplabels,width=150,textvariable=self.grouplabel1_stringvariable,font=font,height=12)
        self.grouplabel2_entry=ctk.CTkEntry(self.row_grouplabels,width=150,textvariable=self.grouplabel2_stringvariable,font=font,height=12)
        self.grouplabel3_entry=ctk.CTkEntry(self.row_grouplabels,width=150,textvariable=self.grouplabel3_stringvariable,font=font,height=12)

        self.colors_buttons=[]
        self.colors_buttons.append(ColorButton(self.row_grouplabels,20,20,color="#27ae60"))
        self.colors_buttons.append(ColorButton(self.row_grouplabels,20,20,color="#e67e22"))
        self.colors_buttons.append(ColorButton(self.row_grouplabels,20,20,color="#808080"))

        if "haplotype_colors" in params:
            for i in range(len(params["haplotype_colors"])):
                if i<3:
                    self.colors_buttons[i].set(params["haplotype_colors"][i])

        self.grouplabel1_entry.grid(row=0,column=1,padx=5,pady=5)
        self.grouplabel2_entry.grid(row=0,column=3,padx=5,pady=5)
        self.update_colors_labels()

        if "group_by" in params: self.update_groupby(params["group_by"])


        # Row4: color by base modification, or SNPs... If base modification
        self.colorby_frame = ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.colorby_frame,10,10,text="Color by: ",font=font).grid(row=0,column=0,pady=5)
        self.color_by="none"
        self.colorby_menu = ctk.CTkOptionMenu(self.colorby_frame,values=["none", "basemod","breakpoints"],font=font,dropdown_font=font,height=12,command=self.update_colorby)
        self.colorby_menu.grid(row=0,column=1,pady=5)
        self.label_unmodified =  ctk.CTkLabel(self.colorby_frame,10,10,text="Unmodified: ",font=font)
        self.colorbutton_unmodified = ColorButton(self.colorby_frame,20,20,color="#0f57e5")
        if "color_unmodified" in params: 
            self.colorbutton_unmodified.set(params["color_unmodified"])
        self.label_basemod1 =  ctk.CTkLabel(self.colorby_frame,10,10,text="Basemod 1: ",font=font)
        self.string_base1 = tkinter.StringVar(self, "C")
        self.entry_base1 = ctk.CTkEntry(self.colorby_frame,30,15,textvariable=self.string_base1,font=font)
        self.string_mod1 = tkinter.StringVar(self, "m")
        self.entry_mod1 = ctk.CTkEntry(self.colorby_frame,30,15,textvariable=self.string_mod1,font=font)
        self.colorbutton_basemod1 = ColorButton(self.colorby_frame,20,20,color="#f40202")
        self.label_basemod2 =  ctk.CTkLabel(self.colorby_frame,10,10,text="Basemod 2: ",font=font)
        self.string_base2 = tkinter.StringVar(self, "C")
        self.entry_base2 = ctk.CTkEntry(self.colorby_frame,30,15,textvariable=self.string_base2,font=font)
        self.string_mod2 = tkinter.StringVar(self, "h")
        self.entry_mod2 = ctk.CTkEntry(self.colorby_frame,30,15,textvariable=self.string_mod2,font=font)
        self.colorbutton_basemod2 = ColorButton(self.colorby_frame,20,20,color="#ffa500")
        self.show_secondbasemod=False
        self.button_add_basemod2 = ctk.CTkButton(self.colorby_frame,20,20,text="Add basemod",command=self.add_basemod,font=font)
        self.button_remove_basemod2 = ctk.CTkButton(self.colorby_frame,20,20,text="Remove basemod",command=self.remove_basemod,font=font)

        if "basemods" in params:
            self.string_base1.set(params["basemods"][0][0])
            self.string_mod1.set(params["basemods"][0][1])
            self.colorbutton_basemod1.set(params["basemods"][0][2])
            self.show_secondbasemod = len(params["basemods"])>1
            if self.show_secondbasemod:
                self.string_base2.set(params["basemods"][1][0])
                self.string_mod2.set(params["basemods"][1][1])
                self.colorbutton_basemod2.set(params["basemods"][1][2])

        # color by breakpoints
        self.breakpoints_entry = EntryPath(self.colorby_frame,label="Breakpoints file: ",entry_width=100)
        if "breakpoints_file" in params: self.breakpoints_entry.set(params["breakpoints_file"])

        if "color_by" in params:  self.update_colorby(params["color_by"])



        self.row1.grid(row=0,column=0,pady=0,sticky="ew")
        self.row2.grid(row=1,column=0,pady=0,sticky="ew")
        self.row3.grid(row=2,column=0,pady=0,sticky="ew")
        self.colorby_frame.grid(row=4,column=0,pady=0,sticky="ew")

    def update_groupby(self,choice):
        if choice!=self.group_by:
            self.group_by=choice
            self.groupby_menu.set(choice)
            if choice=="haplotype":
                self.showunphased_label.grid(row=0,column=2,padx=(5,1),pady=5)
                self.showunphased_checkbox.grid(row=0,column=3,padx=(1,5),pady=5)
                self.exchangehaplotypes_label.grid(row=0,column=4,padx=(5,1),pady=5)
                self.exchangehaplotypes_checkbox.grid(row=0,column=5,padx=(1,5),pady=5)
                self.colors_label.grid(row=0,column=6,padx=(5,1),pady=5)
                self.colors_checkbox.grid(row=0,column=7,padx=(1,5),pady=5)
                self.row_grouplabels.grid(row=3,column=0,pady=0,sticky="ew")
                if self.colors_checkbox.get()==1:
                    if self.showunphased_checkbox.get()==1:
                        self.colors_buttons[2].grid(row=0,column=6,padx=(0,8),pady=5)

            else:
                self.showunphased_label.grid_forget()
                self.showunphased_checkbox.grid_forget()
                self.exchangehaplotypes_label.grid_forget()
                self.exchangehaplotypes_checkbox.grid_forget()
                self.colors_label.grid_forget()
                self.colors_checkbox.grid_forget()
                self.row_grouplabels.grid_forget()
                self.colors_buttons[2].grid_forget()
                #for i in range(3):
                #    self.colors_buttons[i].grid_forget()

    def update_colorby(self,choice):
        if choice!=self.color_by:
            self.color_by=choice
            self.colorby_menu.set(choice)
            if choice=="basemod":
                self.breakpoints_entry.grid_forget()
                self.label_unmodified.grid(row=0,column=2,padx=(5,1),pady=5)
                self.colorbutton_unmodified.grid(row=0,column=3,padx=(1,5),pady=5)
                self.label_basemod1.grid(row=0,column=4,padx=(5,1),pady=5)
                self.entry_base1.grid(row=0,column=5,padx=(1,5),pady=5)
                self.entry_mod1.grid(row=0,column=6,padx=(5,1),pady=5)
                self.colorbutton_basemod1.grid(row=0,column=7,padx=(5,5),pady=5)
                if self.show_secondbasemod:
                    self.button_remove_basemod2.grid(row=1,column=2,padx=(5,1),pady=5,columnspan=2)
                    self.label_basemod2.grid(row=1,column=4,padx=(5,1),pady=5)
                    self.entry_base2.grid(row=1,column=5,padx=(1,5),pady=5)
                    self.entry_mod2.grid(row=1,column=6,padx=(5,1),pady=5)
                    self.colorbutton_basemod2.grid(row=1,column=7,padx=(5,5),pady=5)
                else:
                    self.button_add_basemod2.grid(row=1,column=2,padx=(5,1),pady=5,columnspan=2)
                


            else:
                self.label_unmodified.grid_forget()
                self.colorbutton_unmodified.grid_forget()
                self.label_basemod1.grid_forget()
                self.entry_base1.grid_forget()
                self.entry_mod1.grid_forget()
                self.colorbutton_basemod1.grid_forget()
                self.label_basemod2.grid_forget()
                self.entry_base2.grid_forget()
                self.entry_mod2.grid_forget()
                self.colorbutton_basemod2.grid_forget()
                self.button_add_basemod2.grid_forget()
                self.button_remove_basemod2.grid_forget()
                if choice=="breakpoints":
                    self.breakpoints_entry.grid(row=0,column=2,padx=5,pady=5)

    def add_basemod(self):
        self.show_secondbasemod=True
        self.button_add_basemod2.grid_forget()
        self.button_remove_basemod2.grid(row=1,column=2,padx=(5,1),pady=5,columnspan=2)
        self.label_basemod2.grid(row=1,column=4,padx=(5,1),pady=5)
        self.entry_base2.grid(row=1,column=5,padx=(1,5),pady=5)
        self.entry_mod2.grid(row=1,column=6,padx=(5,1),pady=5)
        self.colorbutton_basemod2.grid(row=1,column=7,padx=(1,5),pady=5)

    def remove_basemod(self):
        self.show_secondbasemod=False
        self.button_add_basemod2.grid(row=1,column=2,padx=(5,1),pady=5,columnspan=2)
        self.button_remove_basemod2.grid_forget()
        self.label_basemod2.grid_forget()
        self.entry_base2.grid_forget()
        self.entry_mod2.grid_forget()
        self.colorbutton_basemod2.grid_forget()

    def update_colors_labels(self):
        if self.colors_checkbox.get()==1:
            self.colors_buttons[0].grid(row=0,column=2,padx=(0,8),pady=5)
            self.colors_buttons[1].grid(row=0,column=4,padx=(0,8),pady=5)
            if self.showunphased_checkbox.get()==1:
                self.colors_buttons[2].grid(row=0,column=6,padx=(0,8),pady=5)
            else:
                self.colors_buttons[2].grid_forget()
        else:
            for i in range(3):
                self.colors_buttons[i].grid_forget()
        if self.showunphased_checkbox.get()==1:
            self.grouplabel3_entry.grid(row=0,column=5,padx=5,pady=5)
        else:
             self.grouplabel3_entry.grid_forget()

    def scalemenu_callback(self,choice):
        if choice!=self.scale:
            if choice=="auto per region": 
                self.scalepos_menu.configure(values=["corner all"])
                self.scalepos_menu.set("corner all")
            elif self.scale=="auto per region": 
                self.scalepos_menu.configure(values=["left","corner","corner all","none"])
                self.scalepos_menu.set("left")
            
            if choice=="custom": 
                self.max_entry.configure(state="normal")
                self.max_stringvariable.set("1")
            else:
                self.max_entry.configure(state="disabled")
                self.max_stringvariable.set("auto")
            self.scale=choice

    def get_params(self):
        params = {}
        params["bam"] =  self.file_entry.get()
        params["label"] = self.label_entry.get()
        params["label_rotate"] = self.label_switch.get()==1
        params["group_by"] = self.group_by
        if params["group_by"]!="none":
            params["exchange_haplotypes"] = self.exchangehaplotypes_checkbox.get()==1
            params["show_unphased"] = self.showunphased_checkbox.get() ==1
            params["color_haplotypes"] = self.colors_checkbox.get()==1
            params["haplotype_colors"] = [b.cget("fg_color") for b in self.colors_buttons]
            params["group_labels"] = [self.grouplabel1_entry.get(),self.grouplabel2_entry.get(),self.grouplabel3_entry.get()] 
        params["color_by"] = self.color_by
        if params["color_by"]=="basemod":
            params["color_unmodified"] = self.colorbutton_unmodified.get()
            params["basemods"] = [(self.entry_base1.get(),self.entry_mod1.get(),self.colorbutton_basemod1.get())]
            if self.show_secondbasemod:  params["basemods"].append((self.entry_base2.get(),self.entry_mod2.get(),self.colorbutton_basemod2.get()))
        elif params["color_by"] == "breakpoints":
            params["breakpoints_file"] = self.breakpoints_entry.get()

        return params

class MetfreqOptionsFrame(ctk.CTkFrame):
    class BamFrame(ctk.CTkFrame):
        def __init__(self, master,params={},**kwargs):
            super().__init__(master.frame_bams,fg_color=colorB,border_width=1)
            self.master2 = master

            # Up/Down buttons
            with resources.as_file(resources.files(figeno.GUI.data) / "up.png") as infile:
                image_up = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(12,12))
            with resources.as_file(resources.files(figeno.GUI.data) / "down.png") as infile:
                image_down = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(12,12))
            self.up_button=ctk.CTkButton(self,width=5,height=5,text="",image=image_up,command=self.up,font=font,fg_color=colorA,hover_color=colorC)
            self.down_button=ctk.CTkButton(self,5,5,text="",image=image_down,command=self.down,font=font,fg_color=colorA,hover_color=colorC)
            with resources.as_file(resources.files(figeno.GUI.data) / "upup.png") as infile:
                image_upup =ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(12,12))
            with resources.as_file(resources.files(figeno.GUI.data) / "downdown.png") as infile:
                image_downdown = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(12,12))
            self.upup_button=ctk.CTkButton(self,width=5,height=5,text="",image=image_upup,command=self.upup,font=font,fg_color=colorA,hover_color=colorC)
            self.downdown_button=ctk.CTkButton(self,5,5,text="",image=image_downdown,command=self.downdown,font=font,fg_color=colorA,hover_color=colorC)
            self.up_button.grid(row=0,column=0,sticky="n",padx=2,pady=2)
            self.down_button.grid(row=1,column=0,sticky="s",padx=2,pady=3)
            self.upup_button.grid(row=0,column=1,sticky="n",padx=2,pady=2)
            self.downdown_button.grid(row=1,column=1,sticky="s",padx=2,pady=3)

            # Main frame
            self.main_frame = ctk.CTkFrame(self,fg_color="transparent")
            # Row1: File, basemod, min coverage, lw, opacity
            self.row1 = ctk.CTkFrame(self.main_frame,fg_color="transparent")
            self.file_entry = EntryPath(self.row1,label="File: ",entry_width=100)
            if "bam" in params: self.file_entry.set(params["bam"])
            self.file_entry.grid(row=0,column=0,padx=(0,5),pady=3)
            self.string_base = tkinter.StringVar(self.row1, "C")
            if "base" in params: self.string_base.set(params["base"])
            self.entry_base = ctk.CTkEntry(self.row1,30,15,textvariable=self.string_base,font=font)
            self.entry_base.grid(row=0,column=3,padx=5,pady=3)
            self.string_mod = tkinter.StringVar(self.row1, "m")
            if "mod" in params: self.string_mod.set(params["mod"])
            self.entry_mod = ctk.CTkEntry(self.row1,30,15,textvariable=self.string_mod,font=font)
            self.entry_mod.grid(row=0,column=4,padx=5,pady=3)

            ctk.CTkLabel(self.row1,10,10,text="min coverage: ",font=font).grid(row=0,column=5,pady=3,padx=(5,1))
            self.mincov_stringvariable = tkinter.StringVar(self.row1, "6")
            if "min_coverage" in params: self.mincov_stringvariable.set(str(params["min_coverage"]))
            self.entry_mincov=ctk.CTkEntry(self.row1,width=40,textvariable=self.mincov_stringvariable,font=font,height=12)
            self.entry_mincov.grid(row=0,column=6,padx=(1,5),pady=3)

            ctk.CTkLabel(self.row1,10,10,text="lw: ",font=font).grid(row=0,column=7,pady=3,padx=(5,1))
            self.lw_stringvariable = tkinter.StringVar(self.row1, "3.0")
            if "linewidth" in params: self.lw_stringvariable.set(str(params["linewidth"]))
            self.entry_lw=ctk.CTkEntry(self.row1,width=40,textvariable=self.lw_stringvariable,font=font,height=12)
            self.entry_lw.grid(row=0,column=8,padx=(1,5),pady=3)
            ctk.CTkLabel(self.row1,10,10,text="opacity: ",font=font).grid(row=0,column=9,pady=3,padx=(5,1))
            self.opacity_stringvariable = tkinter.StringVar(self.row1, "1.0")
            if "alpha" in params: self.opacity_stringvariable.set(str(params["alpha"]))
            self.entry_opacity=ctk.CTkEntry(self.row1,width=40,textvariable=self.opacity_stringvariable,font=font,height=12)
            self.entry_opacity.grid(row=0,column=10,padx=(1,5),pady=3)


            # Row2: split by haplotype and labels/colors
            self.row2 = ctk.CTkFrame(self.main_frame,fg_color="transparent")
            ctk.CTkLabel(self.row2,10,10,text="Split by haplotype: ",font=font).grid(row=0,column=0,pady=3,padx=(5,0))
            self.split_switch = ctk.CTkSwitch(self.row2,10,10,text="",font=font,command=self.update_split_switch)
            if "split" in params and params["split"]: self.split_switch.select()
            self.split_switch.grid(row=0,column=1,padx=(1,5),pady=3)
            ctk.CTkLabel(self.row2,10,10,text="Label 1: ",font=font).grid(row=0,column=2,pady=3,padx=(8,0))
            self.string_label = tkinter.StringVar(self.row2, "")
            if "labels" in params: self.string_label.set(params["labels"][0])
            self.entry_label = ctk.CTkEntry(self.row2,100,15,textvariable=self.string_label,font=font)
            self.entry_label.grid(row=0,column=3,pady=3,padx=5)
            self.colorbutton=ColorButton(self.row2,20,20,color="#1a7242")
            self.colorbutton.grid(row=0,column=4,pady=3,padx=5)

            self.label2=ctk.CTkLabel(self.row2,10,10,text="Label 2: ",font=font)
            self.string_label2 = tkinter.StringVar(self.row2, "")
            if "labels" in params and len(params["labels"])>1: self.string_label2.set(params["labels"][1])
            self.entry_label2 = ctk.CTkEntry(self.row2,100,15,textvariable=self.string_label2,font=font)
            self.colorbutton2=ColorButton(self.row2,20,20,color="#e67e22")
            if "colors" in params:
                self.colorbutton.set(params["colors"][0])
                if len(params["colors"])>1:
                    self.colorbutton2.set(params["colors"][1])
            self.update_split_switch()


            self.row1.grid(row=0,column=0,sticky="ew")
            self.row2.grid(row=1,column=0,sticky="ew")
            self.main_frame.grid(row=0,column=2,rowspan=2,sticky="we",pady=2)


            #image_delete =ctk.CTkImage(Image.open("trash.png").resize((100,100),Image.Resampling.LANCZOS))
            with resources.as_file(resources.files(figeno.GUI.data) / "delete.png") as infile:
                image_delete =ctk.CTkImage(Image.open(infile).resize((100,100),Image.Resampling.LANCZOS),size=(12,15))
            self.delete_button=ctk.CTkButton(self,12,15,text="",image=image_delete,command=self.destroy,font=("arial",2),fg_color="#e74c3c",hover_color=color2hovercolor("#e74c3c"))
            self.delete_button.grid(row=0,column=100,padx=5,pady=5,sticky="e",rowspan=2)
            self.columnconfigure(50,weight=1)
    

        def update_split_switch(self):
            if self.split_switch.get()==1:
                self.label2.grid(row=0,column=5,padx=(8,0))
                self.entry_label2.grid(row=0,column=6,pady=3,padx=5)
                self.colorbutton2.grid(row=0,column=7,pady=3,padx=5)
            else:
                self.label2.grid_forget()
                self.entry_label2.grid_forget()
                self.colorbutton2.grid_forget()
        
        def destroy(self):
            self.master2.remove_bam(self)
            super().destroy()
        def up(self):
            self.master2.move_bamframe_up(self)
        def down(self):
            self.master2.move_bamframe_down(self)
        def upup(self):
            self.master2.move_bamframe_upup(self)
        def downdown(self):
            self.master2.move_bamframe_downdown(self)

        
        def get_params(self):
            params={}
            params["bam"] = self.file_entry.get()
            params["base"] = self.string_base.get()
            params["mod"] = self.string_mod.get()
            params["min_coverage"] = int(float(self.entry_mincov.get()))
            params["linewidth"] = float(self.lw_stringvariable.get())
            params["alpha"] = float(self.opacity_stringvariable.get())
            params["split"] = self.split_switch.get() ==1
            params["labels"] = [self.entry_label.get(),self.entry_label2.get()]
            params["colors"] = [self.colorbutton.get(),self.colorbutton2.get()]
            return params
    class ModbedFrame(ctk.CTkFrame):
        def __init__(self, master,params={},**kwargs):
            super().__init__(master.frame_modbeds,fg_color=colorB,border_width=1)
            self.master2 = master
            ctk.CTkLabel(self,10,10,text="File: ",font=font).grid(row=0,column=0,pady=5,padx=5)
            self.path_stringvariable = tkinter.StringVar(self, "")
            if "modbed" in params: self.path_stringvariable.set(params["modbed"])
            self.entry_path=ctk.CTkEntry(self,width=100,textvariable=self.path_stringvariable,font=font_small,height=12)
            self.entry_path.grid(row=0,column=1,padx=5,pady=5)
            with resources.as_file(resources.files(figeno.GUI.data) / "open.png") as infile:
                image_open = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,13))
            self.button_path = ctk.CTkButton(self,text="",image=image_open,command=self.select_path,font=font,height=12,width=20)
            self.button_path.grid(row=0,column=2,padx=(1,5),pady=5)
            self.string_base = tkinter.StringVar(self, "C")
            if "base" in params: self.string_base.set(params["base"])
            self.entry_base = ctk.CTkEntry(self,30,15,textvariable=self.string_base,font=font)
            self.entry_base.grid(row=0,column=3,padx=5,pady=5)
            self.string_mod = tkinter.StringVar(self, "m")
            if "mod" in params: self.string_mod.set(params["mod"])
            self.entry_mod = ctk.CTkEntry(self,30,15,textvariable=self.string_mod,font=font)
            self.entry_mod.grid(row=0,column=4,padx=5,pady=5)
            ctk.CTkLabel(self,10,10,text="Label: ",font=font).grid(row=0,column=7,pady=5,padx=(8,0))
            self.string_label = tkinter.StringVar(self, "")
            if "label" in params: self.string_label.set(params["label"])
            self.entry_label = ctk.CTkEntry(self,100,15,textvariable=self.string_label,font=font)
            self.entry_label.grid(row=0,column=8,pady=5,padx=5)
            self.colorbutton = ctk.CTkButton(self,20,20,text="",fg_color="#222222",hover_color=color2hovercolor("#222222"),border_width=1)
            self.colorbutton.configure(command=self.select_color(self.colorbutton))
            if "color" in params:
                self.colorbutton.configure(fg_color=params["color"])
                self.colorbutton.configure(hover_color=color2hovercolor(params["color"]))
            self.colorbutton.grid(row=0,column=9,pady=5,padx=5)

            #image_delete =ctk.CTkImage(Image.open("trash.png").resize((100,100),Image.Resampling.LANCZOS))
            with resources.as_file(resources.files(figeno.GUI.data) / "delete.png") as infile:
                image_delete =ctk.CTkImage(Image.open(infile).resize((100,100),Image.Resampling.LANCZOS),size=(12,15))
            self.delete_button=ctk.CTkButton(self,12,15,text="",image=image_delete,command=self.destroy,font=("arial",2),fg_color="#e74c3c",hover_color=color2hovercolor("#e74c3c"))
            self.delete_button.grid(row=0,column=100,padx=5,pady=5,sticky="e")
            self.columnconfigure(50,weight=1)
        

        def select_path(self):
            result = tkinter.filedialog.askopenfilename()
            if len(result)>0:
                self.path_stringvariable.set(result)
                #self.entry_path.focus_set()
                self.entry_path.xview(tkinter.END)

        def select_color(self,button):
            def select_color_b():
                colorselector = ColorSelector(self,button,default_color=button.cget("fg_color"))
            return select_color_b
        
        def destroy(self):
            self.master2.remove_modbed(self)
            super().destroy()

        
        def get_params(self):
            params={}
            params["modbed"] = self.path_stringvariable.get()
            params["base"] = self.string_base.get()
            params["mod"] = self.string_mod.get()
            params["split"] = self.split_switch.get() ==1
            params["label"] = self.entry_label.get()
            params["color"] = self.colorbutton.cget("fg_color")
            return params

    def __init__(self, master,params={},**kwargs):
        super().__init__(master,fg_color="transparent",border_width=0)

        # Row1: label
        # Frame bams (potentially split by haplotype). Base, mod, label, color, linewidth
        # Frame modbeds

        # Row1: label
        self.row1=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row1,10,10,text="Label: ",font=font).grid(row=0,column=0,pady=5)
        self.label_stringvariable = tkinter.StringVar(self, "Methylation freq")
        if "label" in params: self.label_stringvariable.set(params["label"])
        self.label_entry=ctk.CTkEntry(self.row1,width=200,textvariable=self.label_stringvariable,font=font,height=12)
        self.label_entry.grid(row=0,column=1,padx=5,pady=5)
        ctk.CTkLabel(self.row1,10,10,text="Rotate label:",font=font).grid(row=0,column=2,padx=5,pady=5)
        self.label_switch = ctk.CTkSwitch(self.row1,10,10,text="",font=font)
        if (not "label_rotate" in params) or params["label_rotate"]: self.label_switch.select()
        self.label_switch.grid(row=0,column=3,padx=5,pady=5)

        # Frame bams
        self.frame_bams = ctk.CTkFrame(self,corner_radius=5,border_width=1,fg_color=colorC)
        self.frame_bams.columnconfigure(2,weight=1)
        ctk.CTkLabel(self.frame_bams,10,10,text="Bam files",font=font).grid(row=0,column=0,pady=5,padx=5)
        self.addbam_button = ctk.CTkButton(self.frame_bams,20,20,text="Add bam",command=self.add_bam)
        self.addbam_button.grid(row=0,column=1,pady=5,padx=5)
        self.bam_frames=[]
        if "bams" in params:
            for p in params["bams"]:
                self.bam_frames.append(self.BamFrame(self,p))
            self.update_bams()

        # Frame modbeds
        self.frame_modbeds = ctk.CTkFrame(self,corner_radius=5,border_width=1,fg_color=colorC)
        self.frame_modbeds.columnconfigure(2,weight=1)
        ctk.CTkLabel(self.frame_modbeds,10,10,text="Modbed files",font=font).grid(row=0,column=0,pady=5,padx=5)
        self.addmodbed_button = ctk.CTkButton(self.frame_modbeds,20,20,text="Add modbed",command=self.add_modbed)
        self.addmodbed_button.grid(row=0,column=1,pady=5,padx=5)
        self.modbed_frames=[]
        if "modbeds" in params:
            for p in params["modbeds"]:
                self.modbed_frames.append(self.ModbedFrame(self,p))
                



        
        self.row1.grid(row=0,column=0,pady=0,sticky="ew")
        self.frame_bams.grid(row=1,column=0,pady=0,sticky="ew")
        self.frame_modbeds.grid(row=2,column=0,pady=2,sticky="ew")
        self.grid_columnconfigure(0,weight=1)

        

    def add_bam(self):
        self.bam_frames.append(self.BamFrame(self))
        self.update_bams()

    def remove_bam(self,bamframe):
        self.bam_frames.remove(bamframe)

    def move_bamframe_up(self,frame):
        for i in range(1,len(self.bam_frames)):
            if self.bam_frames[i]==frame:
               self.bam_frames[i-1],self.bam_frames[i] = self.bam_frames[i],self.bam_frames[i-1]
               break
        self.update_bams()

    def move_bamframe_down(self,frame):
        for i in range(len(self.bam_frames)-1):
            if self.bam_frames[i]==frame:
                self.bam_frames[i+1],self.bam_frames[i] = self.bam_frames[i],self.bam_frames[i+1]
                break
        self.update_bams()
    def move_bamframe_upup(self,frame):
        self.bam_frames.remove(frame)
        self.bam_frames = [frame] + self.bam_frames
        self.update_bams()

    def move_bamframe_downdown(self,frame):
        self.bam_frames.remove(frame)
        self.bam_frames =  self.bam_frames + [frame]
        self.update_bams()

    def update_bams(self):
        for i in range(len(self.bam_frames)):
            self.bam_frames[i].grid(row=1+i,column=0,columnspan=100,pady=2,padx=3,sticky="ew")
    def add_modbed(self):
        self.modbed_frames.append(self.ModbedFrame(self))
        self.update_modbeds()
    def remove_modbed(self,modbedframe):
        self.modbed_frames.remove(modbedframe)

    def update_modbeds(self):
        for i in range(len(self.modbed_frames)):
            self.modbed_frames[i].grid(row=1+i,column=0,columnspan=100,pady=2,padx=3,sticky="ew")

    def get_params(self):
        params = {}
        params["label"] = self.label_entry.get()
        params["label_rotate"] = self.label_switch.get()==1
        params["bams"] = []
        for x in self.bam_frames: params["bams"].append(x.get_params())
        params["modbeds"] = []
        for x in self.modbed_frames: params["modbeds"].append(x.get_params())

        return params


class HiCOptionsFrame(ctk.CTkFrame):
    def __init__(self, master,params={},**kwargs):
        super().__init__(master,fg_color="transparent",border_width=0)


        # Row1: file path
        self.row1=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        self.entry_file = EntryPath(self.row1,label="File: ",entry_width=300)
        if "file" in params: self.entry_file.set(params["file"])
        self.entry_file.grid(row=0,column=0,padx=(0,5),pady=5)

        # Row2: label
        self.row2=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row2,10,10,text="Label: ",font=font).grid(row=0,column=0,pady=5)
        self.label_stringvariable = tkinter.StringVar(self, "")
        if "label" in params: self.label_stringvariable.set(params["label"])
        self.label_entry=ctk.CTkEntry(self.row2,width=200,textvariable=self.label_stringvariable,font=font,height=12)
        self.label_entry.grid(row=0,column=1,padx=5,pady=5)
        ctk.CTkLabel(self.row2,10,10,text="Rotate label:",font=font).grid(row=0,column=2,padx=5,pady=5)
        self.label_switch = ctk.CTkSwitch(self.row2,10,10,text="",font=font)
        if "label_rotate" in params and params["label_rotate"]: self.label_switch.select()
        self.label_switch.grid(row=0,column=3,padx=5,pady=5)

        # Row3: colormap, pixel border, upside down
        self.row3=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row3,10,10,text="Color map: ",font=font).grid(row=0,column=0,pady=5)
        self.colormap_menu = ctk.CTkOptionMenu(self.row3,values=["Red", "Heat"],font=font,dropdown_font=font,height=12)
        if "color_map" in params: self.colormap_menu.set(params["color_map"])
        self.colormap_menu.grid(row=0,column=1,pady=5)
        ctk.CTkLabel(self.row3,10,10,text="Pixel border: ",font=font).grid(row=0,column=2,pady=5,padx=(10,0))
        self.pixelborder_switch = ctk.CTkSwitch(self.row3,10,10,text="",font=font)
        self.pixelborder_switch.grid(row=0,column=3,pady=5,padx=(0,5))
        if ("pixel_border" in params) and params["pixel_border"]: self.pixelborder_switch.select()
        ctk.CTkLabel(self.row3,10,10,text="Upside down: ",font=font).grid(row=0,column=8,pady=5,padx=(10,0))
        self.upsidedown_switch = ctk.CTkSwitch(self.row3,10,10,text="",font=font)
        self.upsidedown_switch.grid(row=0,column=9,pady=5,padx=(0,5))
        if "upside_down" in params and params["upside_down"]: self.upsidedown_switch.select()
        

        # Row4: maxdist, across reg
        self.row4=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row4,10,10,text="Max dist (kb): ",font=font).grid(row=0,column=0,pady=5)
        self.maxdist_stringvariable = tkinter.StringVar(self, "700")
        if "max_dist" in params: self.maxdist_stringvariable.set(str(params["max_dist"]))
        self.entry_maxdist=ctk.CTkEntry(self.row4,width=60,textvariable=self.maxdist_stringvariable,font=font,height=12)
        self.entry_maxdist.grid(row=0,column=1,padx=5,pady=5)
        ctk.CTkLabel(self.row4,10,10,text="Extend: ",font=font).grid(row=0,column=2,pady=5,padx=(10,0))
        self.extend_switch = ctk.CTkSwitch(self.row4,10,10,text="",font=font)
        self.extend_switch.grid(row=0,column=3,pady=5,padx=(0,5))
        if (not "extend" in params) or params["extend"]: self.extend_switch.select()
        ctk.CTkLabel(self.row4,10,10,text="Interactions across regions: ",font=font).grid(row=0,column=4,pady=5,padx=(10,0))
        self.acrossreg_switch = ctk.CTkSwitch(self.row4,10,10,text="",font=font,command=self.update_across_regions)
        self.acrossreg_switch.grid(row=0,column=5,pady=5,padx=(0,5))
        if (not "interactions_across_regions" in params) or params["interactions_across_regions"]: self.acrossreg_switch.select()
        self.double_label=ctk.CTkLabel(self.row4,10,10,text="Double across regions: ",font=font)
        self.double_switch = ctk.CTkSwitch(self.row4,10,10,text="",font=font)
        if (not "double_interactions_across_regions" in params) or params["double_interactions_across_regions"]: self.double_switch.select()
        self.update_across_regions()
        
        

        self.row1.grid(row=0,column=0,pady=0,sticky="ew")
        self.row2.grid(row=1,column=0,pady=0,sticky="ew")
        self.row3.grid(row=2,column=0,pady=0,sticky="ew")
        self.row4.grid(row=3,column=0,pady=0,sticky="ew")
        self.grid_columnconfigure(0,weight=1)

    def update_across_regions(self):
        if self.acrossreg_switch.get()==1:
            self.double_label.grid(row=0,column=6,pady=5,padx=(10,0))
            self.double_switch.grid(row=0,column=7,pady=5,padx=(2,0))
        else:
            self.double_label.grid_forget()
            self.double_switch.grid_forget()


    def get_params(self):
        params = {}
        params["file"] =  self.entry_file.get()
        params["label"] = self.label_entry.get()
        params["label_rotate"] = self.label_switch.get()==1
        params["color_map"] = self.colormap_menu.get()
        params["pixel_border"] = self.pixelborder_switch.get()==1
        params["max_dist"] = int(self.entry_maxdist.get())
        params["interactions_across_regions"] = self.acrossreg_switch.get()==1
        if self.acrossreg_switch.get()==1:
            params["double_interactions_across_regions"] = self.double_switch.get()==1
        params["extend"] = self.extend_switch.get()==1
        params["upside_down"] = self.upsidedown_switch.get()==1

        return params
    
class SVOptionsFrame(ctk.CTkFrame):
    def __init__(self, master,params={},**kwargs):
        super().__init__(master,fg_color="transparent",border_width=0)

        #TODO: allow user-defined SVs instead of vcf ??
        # Row1: vcf
        self.row1=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        self.file_entry = EntryPath(self.row1,"File: ",entry_width=300)
        if "vcf" in params: self.file_entry.set(params["vcf"])
        self.file_entry.grid(row=0,column=1,pady=5,padx=(0,5))

        # Row2: label
        self.row2=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row2,10,10,text="Label: ",font=font).grid(row=0,column=0,pady=5)
        self.label_stringvariable = tkinter.StringVar(self, "BP")
        if "label" in params: self.label_stringvariable.set(params["label"])
        self.label_entry=ctk.CTkEntry(self.row2,width=200,textvariable=self.label_stringvariable,font=font,height=12)
        self.label_entry.grid(row=0,column=1,padx=5,pady=5)
        ctk.CTkLabel(self.row2,10,10,text="Rotate label:",font=font).grid(row=0,column=2,padx=5,pady=5)
        self.label_switch = ctk.CTkSwitch(self.row2,10,10,text="",font=font)
        if (not "label_rotate" in params) or params["label_rotate"]: self.label_switch.select()
        self.label_switch.grid(row=0,column=3,padx=5,pady=5)


        
        # Row 3: Colors
        self.row3=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row3,10,10,text="DEL-like",font=font).grid(row=0,column=0,padx=(0,2),pady=5)
        self.buttoncolor_del = ColorButton(self.row3,20,20,color="#4a69bd")
        if "color_del" in params: self.buttoncolor_del.set(params["color_del"])
        self.buttoncolor_del.grid(row=0,column=1,padx=2,pady=5)
        ctk.CTkLabel(self.row3,10,10,text="DUP-like",font=font).grid(row=0,column=2,padx=(8,2),pady=5)
        self.buttoncolor_dup = ColorButton(self.row3,20,20,color="#e55039")
        if "color_dup" in params: self.buttoncolor_dup.set(params["color_dup"])
        self.buttoncolor_dup.grid(row=0,column=3,padx=2,pady=5)
        ctk.CTkLabel(self.row3,10,10,text="T2T",font=font).grid(row=0,column=4,padx=(8,2),pady=5)
        self.buttoncolor_t2t = ColorButton(self.row3,20,20,color="#8e44ad")
        if "color_t2t" in params: self.buttoncolor_t2t.set(params["color_h2h"])
        self.buttoncolor_t2t.grid(row=0,column=5,padx=2,pady=5)
        ctk.CTkLabel(self.row3,10,10,text="H2H",font=font).grid(row=0,column=6,padx=(8,2),pady=5)
        self.buttoncolor_h2h = ColorButton(self.row3,20,20,color="#8e44ad")
        if "color_h2h" in params: self.buttoncolor_h2h.set(params["color_t2t"])
        self.buttoncolor_h2h.grid(row=0,column=7,padx=2,pady=5)
        ctk.CTkLabel(self.row3,10,10,text="Translocation",font=font).grid(row=0,column=8,padx=(8,2),pady=5)
        self.buttoncolor_trans = ColorButton(self.row3,20,20,color="#27ae60")
        if "color_trans" in params: self.buttoncolor_trans.set(params["color_trans"])
        self.buttoncolor_trans.grid(row=0,column=9,padx=2,pady=5)
        #self.colorRow = colorFrame(self,"Color")        
        
        self.row1.grid(row=0,column=0,pady=0,sticky="ew")
        self.row2.grid(row=1,column=0,pady=0,sticky="ew")
        self.row3.grid(row=2,column=0,pady=0,sticky="ew")

    def get_params(self):
        params =  {"vcf":self.file_entry.get()}
        params["label"] = self.label_entry.get()
        params["label_rotate"] = self.label_switch.get() == 1
        params["color_del"] = self.buttoncolor_del.get()
        params["color_dup"] = self.buttoncolor_dup.get()
        params["color_h2h"] = self.buttoncolor_h2h.get()
        params["color_t2t"] = self.buttoncolor_t2t.get()
        params["color_trans"] = self.buttoncolor_trans.get()

        return params
    


class CopynumberOptionsFrame(ctk.CTkFrame):
    def __init__(self, master,params={},**kwargs):
        super().__init__(master,fg_color="transparent",border_width=0)

        # Row1: Ratios and CNA file
        self.row1=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        self.ratio_entry = EntryPath(self.row1,"Ratios file: ",entry_width=140)
        if "freec_ratios" in params: self.ratio_entry.set(params["freec_ratios"])
        self.ratio_entry.grid(row=0,column=0,pady=5,padx=(0,5))

        self.cnas_entry = EntryPath(self.row1,"CNAs file: ",entry_width=140)
        if "freec_CNAs" in params: self.cnas_entry.set(params["freec_CNAs"])
        self.cnas_entry.grid(row=0,column=1,pady=5,padx=(3,5))

        self.purple_entry = EntryPath(self.row1,"Purple CN file: ",entry_width=140)
        if "purple_cn" in params: self.purple_entry.set(params["purple_cn"])
        self.purple_entry.grid(row=0,column=2,pady=5,padx=(3,5))

        # Row2: label
        self.row2=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row2,10,10,text="Label: ",font=font).grid(row=0,column=0,pady=5)
        self.label_stringvariable = tkinter.StringVar(self, "CN")
        if "label" in params: self.label_stringvariable.set(params["label"])
        self.label_entry=ctk.CTkEntry(self.row2,width=200,textvariable=self.label_stringvariable,font=font,height=12)
        self.label_entry.grid(row=0,column=1,padx=5,pady=5)
        ctk.CTkLabel(self.row2,10,10,text="Rotate label:",font=font).grid(row=0,column=2,padx=5,pady=5)
        self.label_switch = ctk.CTkSwitch(self.row2,10,10,text="",font=font)
        if (not "label_rotate" in params) or params["label_rotate"]: self.label_switch.select()
        self.label_switch.grid(row=0,column=3,padx=5,pady=5)
        ctk.CTkLabel(self.row2,10,10,text="Genes:",font=font).grid(row=0,column=4,padx=(10,1),pady=5)
        self.genes_stringvariable = tkinter.StringVar(self, "")
        if "genes_highlighted" in params: self.genes_stringvariable.set(",".join(params["genes_highlighted"]))
        self.genes_entry = ctk.CTkEntry(self.row2,width=200,textvariable=self.genes_stringvariable,font=font,height=12)
        self.genes_entry.grid(row=0,column=5,padx=5,pady=5)

        
        # Row 3: Min/Max and Colors and grid
        self.row3=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row3,10,10,text="Min CN: ",font=font).grid(row=0,column=0,padx=(0,2),pady=5)
        self.mincn_stringvariable = tkinter.StringVar(self, "")
        if "min_cn" in params: self.mincn_stringvariable.set(str(params["min_cn"]))
        self.mincn_entry = ctk.CTkEntry(self.row3,width=60,font=font,textvariable=self.mincn_stringvariable)
        self.mincn_entry.grid(row=0,column=1,padx=(1,5),pady=5)

        ctk.CTkLabel(self.row3,10,10,text="Max CN: ",font=font).grid(row=0,column=2,padx=(0,2),pady=5)
        self.maxcn_stringvariable = tkinter.StringVar(self, "")
        if "max_cn" in params: self.maxcn_stringvariable.set(str(params["max_cn"]))
        self.maxcn_entry = ctk.CTkEntry(self.row3,width=60,font=font,textvariable=self.maxcn_stringvariable)
        self.maxcn_entry.grid(row=0,column=3,padx=(1,5),pady=5)


        ctk.CTkLabel(self.row3,10,10,text="Normal",font=font).grid(row=0,column=4,padx=(8,2),pady=5)
        self.buttoncolor_normal = ColorButton(self.row3,20,20,color="#000000")
        if "color_normal" in params: self.buttoncolor_normal.set(params["color_normal"])
        self.buttoncolor_normal.grid(row=0,column=5,padx=2,pady=5)

        ctk.CTkLabel(self.row3,10,10,text="Loss",font=font).grid(row=0,column=6,padx=(8,2),pady=5)
        self.buttoncolor_loss = ColorButton(self.row3,20,20,color="#4a69bd")
        if "color_loss" in params: self.buttoncolor_loss.set(params["color_loss"])
        self.buttoncolor_loss.grid(row=0,column=7,padx=2,pady=5)

        ctk.CTkLabel(self.row3,10,10,text="Gain",font=font).grid(row=0,column=8,padx=(8,2),pady=5)
        self.buttoncolor_gain = ColorButton(self.row3,20,20,color="#e55039")
        if "color_gain" in params: self.buttoncolor_gain.set(params["color_gain"])
        self.buttoncolor_gain.grid(row=0,column=9,padx=2,pady=5)

        ctk.CTkLabel(self.row3,10,10,text="Grid: ",font=font).grid(row=0,column=10,padx=(10,2),pady=5)
        self.grid_switch = ctk.CTkSwitch(self.row3,10,10,text="")

        if (not "grid" in params) or (params["grid"]): 
            self.grid_switch.select()
        self.grid_switch.grid(row=0,column=11,padx=2,pady=5)
        
        
        self.row1.grid(row=0,column=0,pady=0,sticky="ew")
        self.row2.grid(row=1,column=0,pady=0,sticky="ew")
        self.row3.grid(row=2,column=0,pady=0,sticky="ew")

    def get_params(self):
        params =  {}
        if self.ratio_entry.get()!="": params["freec_ratios"] = self.ratio_entry.get()
        if self.cnas_entry.get()!="": params["freec_CNAs"] = self.cnas_entry.get()
        if self.purple_entry.get()!="": params["purple_cn"] = self.purple_entry.get()
        if self.mincn_stringvariable.get()!="": params["min_cn"] = float(self.mincn_stringvariable.get())
        if self.maxcn_stringvariable.get()!="": params["max_cn"] = float(self.maxcn_stringvariable.get())
        if self.genes_stringvariable.get()!="": params["genes_highlighted"] = self.genes_stringvariable.get().split(",")
        params["color_normal"] = self.buttoncolor_normal.get()
        params["color_loss"] = self.buttoncolor_loss.get()
        params["color_gain"] = self.buttoncolor_gain.get()
        params["grid"] = self.grid_switch.get()==1
        params["label"] = self.label_entry.get()
        params["label_rotate"] = self.label_switch.get() == 1

        return params



class AseOptionsFrame(ctk.CTkFrame):
    def __init__(self, master,params={},**kwargs):
        super().__init__(master,fg_color="transparent",border_width=0)

        # Row1: RNA and DNA files
        self.row1=ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        self.rna_entry = EntryPath(self.row1,label="RNA file: ",entry_width=150)
        if "ase_file" in params: self.rna_entry.set(params["ase_file"])
        self.rna_entry.grid(column=0,row=0,pady=5)
        self.dna_entry = EntryPath(self.row1,label="DNA file: ",entry_width=150)
        if "vcf_DNA" in params: self.dna_entry.set(params["vcf_DNA"])
        self.dna_entry.grid(column=1,row=0,pady=5,padx=8)

        # Row2: Other params
        self.row2 = ctk.CTkFrame(self,corner_radius=0,border_width=0,fg_color="transparent")
        ctk.CTkLabel(self.row2,10,10,text="Min depth: ",font=font).grid(column=0,row=0,pady=5)
        self.mindepth_stringvar = tkinter.StringVar(self, "6")
        if "min_depth" in params: self.mindepth_stringvar.set(str(params["min_depth"]))
        self.mindepth_entry = ctk.CTkEntry(self.row2,50,10,textvariable=self.mindepth_stringvar,font=font)
        self.mindepth_entry.grid(row=0,column=1,pady=5,padx=(1,5))
        ctk.CTkLabel(self.row2,10,10,text="Only exonic: ",font=font).grid(column=2,row=0,pady=5,padx=5)
        self.onlyexonic_checkbox = ctk.CTkCheckBox(self.row2,width=20,height=20,checkbox_width=20,checkbox_height=20,text="")
        if "only_exonic" in params and params["only_exonic"]: self.onlyexonic_checkbox.select()
        self.onlyexonic_checkbox.grid(column=3,row=0,pady=5,padx=(0,5))
        ctk.CTkLabel(self.row2,10,10,text="Colors: ",font=font).grid(column=4,row=0,pady=5,padx=(5,1))
        self.colorbutton1=ColorButton(self.row2,color="#e55039")
        if "color1" in params: self.colorbutton1.set(params["color1"])
        self.colorbutton1.grid(column=5,row=0,pady=5,padx=(0,1))
        self.colorbutton2=ColorButton(self.row2,color="#4a69bd")
        if "color2" in params: self.colorbutton2.set(params["color2"])
        self.colorbutton2.grid(column=6,row=0,pady=5,padx=5)
        ctk.CTkLabel(self.row2,10,10,text="Grid: ",font=font).grid(column=7,row=0,pady=5,padx=(8,0))
        self.grid_checkbox = ctk.CTkCheckBox(self.row2,width=20,height=20,checkbox_width=20,checkbox_height=20,text="")
        if "grid" in params and params["grid"]: self.grid_checkbox.select()
        self.grid_checkbox.grid(column=8,row=0,pady=5,padx=5)

        self.row1.grid(row=0,column=0,pady=0,sticky="ew")
        self.row2.grid(row=1,column=0,pady=0,sticky="ew")


    def get_params(self):
        params = {}
        params["ase_file"] = self.rna_entry.get()
        params["vcf_DNA"] = self.dna_entry.get()
        params["min_depth"] = int(self.mindepth_entry.get())
        params["only_exonic"] = self.onlyexonic_checkbox.get()==1
        params["color1"] = self.colorbutton1.get()
        params["color2"] = self.colorbutton2.get()
        params["grid"] = self.grid_checkbox.get()==1
        return params