import customtkinter as ctk
import tkinter
from PIL import Image
import importlib_resources as resources
from figeno.GUI.gui_colors import ColorButton, color2hovercolor
import figeno.GUI.data
font_title=("inter",16,"bold")
font=("inter",16)
font_small=("inter",10)



colorE="gray35"
colorD="gray25"
colorC="gray16"
colorB="gray13"
colorA="gray10"

class RegionsFrame(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master,fg_color=colorA,corner_radius=3,border_width=1)
        self.region_frames=[]
        self.unused_region_frames=[]

        self.label = ctk.CTkLabel(self,text="Regions",font=font_title)
        self.label.grid(row=0,column=0,padx=10,pady=5,sticky="w")
        with resources.as_file(resources.files(figeno.GUI.data) / "plus.png") as infile:
            image_plus = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(17,17))
        self.button =ctk.CTkButton(self, text="Add region", command=self.add_region_frame,font=font,image=image_plus,height=16)
        self.button.grid(row=0, column=1, padx=10, pady=(6,1), sticky="w")
        self.button_allchr =ctk.CTkButton(self, text="Add all chromosomes", command=self.add_all_chromosomes,font=font,image=image_plus,height=16)
        self.button_allchr.grid(row=0, column=2, padx=10, pady=(6,1), sticky="w")
        #self.button2 =ctk.CTkButton(self, text="Display regions")
        #self.button2.grid(row=0, column=2, padx=10, pady=0, sticky="w")

        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=0)
        self.columnconfigure(2, weight=1)
        self.add_region_frame()

    def add_region_frame(self,params={},update_grid=True):
        if len(self.unused_region_frames)>0:
            region_frame = self.unused_region_frames.pop()
            region_frame.reload(params)
        else:
            region_frame = RegionFrame(self,params)
        self.region_frames.append(region_frame)
        if update_grid:
            self.update_regionframes_grid()
    def add_all_chromosomes(self):
        colors=["#98671F","#65661B","#969833","#CE151D","#FF1A25","#FF0BC8","#FFCBCC","#FF9931","#FFCC3A","#FCFF44","#C4FF40","#00FF3B",
            "#2F7F1E","#2800C6","#6A96FA","#98CAFC","#00FEFD","#C9FFFE","#9D00C6","#D232FA","#956DB5","#5D5D5D","#989898","#CBCBCB"]
        chromosomes = [str(x) for x in range(1,23)] + ["X","Y"]
        for i in range(len(chromosomes)):
            self.add_region_frame({"chr":chromosomes[i],"color":colors[i]},update_grid=False)
        self.update_regionframes_grid()

    def remove_region_frame(self,regionframe):
        regionframe.grid_forget()
        self.unused_region_frames.append(regionframe)
        self.region_frames.remove(regionframe)

    def move_frame_up(self,frame):
        for i in range(1,len(self.region_frames)):
            if self.region_frames[i]==frame:
                self.region_frames[i-1],self.region_frames[i] = self.region_frames[i],self.region_frames[i-1]
                break
        self.update_regionframes_grid()

    def move_frame_down(self,frame):
        for i in range(len(self.region_frames)-1):
            if self.region_frames[i]==frame:
                self.region_frames[i+1],self.region_frames[i] = self.region_frames[i],self.region_frames[i+1]
                break
        self.update_regionframes_grid()

    def upup(self,frame):
        self.region_frames.remove(frame)
        self.region_frames = [frame] + self.region_frames
        self.update_regionframes_grid()
    def downdown(self,frame):
        self.region_frames.remove(frame)
        self.region_frames = self.region_frames + [frame]
        self.update_regionframes_grid()

    def update_regionframes_grid(self):
        for i in range(len(self.region_frames)):
            pady=(2,2)
            if i==0: pady =  (4,pady[1])
            if i==len(self.region_frames)-1: pady =  (pady[0],5)
            self.region_frames[i].grid(row=i+1, column=0, padx=10, pady=pady, sticky="new",columnspan=3)

    def reload(self,params_list=[]):
        if len(params_list)==0:
            params_list.append({})
        n_previous_regions = len(self.region_frames)
        n_new_regions = len(params_list)
        if n_previous_regions>n_new_regions:
            for i in range(n_previous_regions-n_new_regions):
                region_frame = self.region_frames.pop()
                region_frame.grid_forget()
                self.unused_region_frames.append(region_frame)
        for i in range(min(n_new_regions,n_previous_regions)):
            self.region_frames[i].reload(params_list[i])
        if n_new_regions>n_previous_regions:
            for i in range(n_previous_regions,n_new_regions):
                self.add_region_frame(params_list[i],update_grid=False)
        self.update_regionframes_grid()


    def get_regions(self):
        l=[]
        for x in self.region_frames:
           l.append(x.get_region())
        return l

class RegionFrame(ctk.CTkFrame):
    def __init__(self, master,params={}):
        super().__init__(master,fg_color=colorB,border_width=2)

        with resources.as_file(resources.files(figeno.GUI.data) / "up.png") as infile:
            image_up = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,15))
        with resources.as_file(resources.files(figeno.GUI.data) / "down.png") as infile:
            image_down = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,15))
        self.up_button=ctk.CTkButton(self,width=5,height=5,text="",image=image_up,command=self.up,font=font,fg_color=colorA,hover_color=colorC)
        self.down_button=ctk.CTkButton(self,5,5,text="",image=image_down,command=self.down,font=font,fg_color=colorA,hover_color=colorC)
        self.up_button.grid(row=0,column=0,padx=2,pady=2)
        self.down_button.grid(row=1,column=0,padx=2,pady=2)

        with resources.as_file(resources.files(figeno.GUI.data) / "upup.png") as infile:
            image_upup =ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,15))
        with resources.as_file(resources.files(figeno.GUI.data) / "downdown.png") as infile:
            image_downdown = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,15))
        self.upup_button=ctk.CTkButton(self,width=5,height=5,text="",image=image_upup,command=self.upup,font=font,fg_color=colorA,hover_color=colorC)
        self.downdown_button=ctk.CTkButton(self,5,5,text="",image=image_downdown,command=self.downdown,font=font,fg_color=colorA,hover_color=colorC)
        self.upup_button.grid(row=0,column=1,padx=2,pady=2)
        self.downdown_button.grid(row=1,column=1,padx=2,pady=2)

        self.chr_label=ctk.CTkLabel(self,10,20,text="chr:",font=font)
        self.chr_label.grid(row=0,column=2,rowspan=2,padx=(10,0))
        self.chr_stringvar = tkinter.StringVar(self, "")
        if "chr" in params: self.chr_stringvar.set(params["chr"])
        self.chr_entry=ctk.CTkEntry(self,50,20,font=font,textvariable=self.chr_stringvar)
        self.chr_entry.grid(row=0,column=3,rowspan=2,padx=(0,10))

        self.start_label=ctk.CTkLabel(self,10,20,text="start:",font=font)
        self.start_label.grid(row=0,column=4,rowspan=2,padx=(10,0))
        self.start_stringvar = tkinter.StringVar(self, "")
        if "start" in params: self.start_stringvar.set(params["start"])
        self.start_entry=ctk.CTkEntry(self,200,20,font=font,textvariable=self.start_stringvar)
        self.start_entry.grid(row=0,column=5,rowspan=2,padx=(0,10))

        self.end_label=ctk.CTkLabel(self,10,10,text="end:",font=font)
        self.end_label.grid(row=0,column=6,rowspan=2,padx=(10,0))
        self.end_stringvar = tkinter.StringVar(self, "")
        if "end" in params: self.end_stringvar.set(str(params["end"]))
        self.end_entry=ctk.CTkEntry(self,200,20,font=font,textvariable=self.end_stringvar)
        self.end_entry.grid(row=0,column=7,rowspan=2,padx=(0,10))

        self.colorbutton=ColorButton(self,25,25,color="#f4a460")
        if "color" in params:
            self.colorbutton.set(params["color"])
        self.colorbutton.grid(row=0,column=8,padx=(15,2),rowspan=2)
       
        with resources.as_file(resources.files(figeno.GUI.data) / "delete.png") as infile:
            image_delete =ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS))
        self.delete_button=ctk.CTkButton(self,10,10,text="",image=image_delete,command=self.removeframe,font=font,fg_color="#e74c3c",hover_color=color2hovercolor("#e74c3c"))
        self.delete_button.grid(row=0,column=9,rowspan=2,padx=10,sticky="e")
        self.grid_columnconfigure(9,weight=1)

    def removeframe(self):
        self.master.remove_region_frame(self)
    def up(self):
        self.master.move_frame_up(self)
    def down(self):
        self.master.move_frame_down(self)
    def upup(self):
        self.master.upup(self)
    def downdown(self):
        self.master.downdown(self)
    
    def reload(self,params={}):
        if "chr" in params: self.chr_stringvar.set(params["chr"])
        else: self.chr_stringvar.set("")
        if "start" in params: self.start_stringvar.set(str(params["start"]))
        else: self.start_stringvar.set("")
        if "end" in params: self.end_stringvar.set(str(params["end"]))
        else: self.end_stringvar.set("")
        if "color" in params: self.colorbutton.set(params["color"])
        else: self.colorbutton.set("#f4a460")

    def get_region(self):
        params =  {"chr":self.chr_entry.get()}
        try: params["start"]=int(self.start_entry.get())
        except: pass
        try: params["end"]=int(self.end_entry.get())
        except: pass
        params["color"]=self.colorbutton.get()
        return params
    


class HighlightsFrame(ctk.CTkFrame):
    def __init__(self, master):
        super().__init__(master,fg_color=colorA,corner_radius=3,border_width=1)
        self.region_frames=[]
        self.unused_region_frames=[]

        self.label = ctk.CTkLabel(self,text="Highlights",font=font_title)
        self.label.grid(row=0,column=0,padx=10,pady=5,sticky="w")
        with resources.as_file(resources.files(figeno.GUI.data) / "plus.png") as infile:
            image_plus = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(17,17))
        self.button =ctk.CTkButton(self, text="Add Highlight", command=self.add_region_frame,font=font,image=image_plus,height=16)
        self.button.grid(row=0, column=1, padx=10, pady=0, sticky="w")
        #self.button2 =ctk.CTkButton(self, text="Display regions")
        #self.button2.grid(row=0, column=2, padx=10, pady=0, sticky="w")

        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=0)
        self.columnconfigure(2, weight=1)

    def add_region_frame(self,params={},update_grid=True):
        if len(self.unused_region_frames)>0:
            region_frame = self.unused_region_frames.pop()
            region_frame.reload(params)
        else:
            region_frame = HighlightFrame(self,params)
        self.region_frames.append(region_frame)
        if update_grid:
            self.update_regionframes_grid()

    def remove_region_frame(self,regionframe):
        regionframe.grid_forget()
        self.unused_region_frames.append(regionframe)
        self.region_frames.remove(regionframe)

    def move_frame_up(self,frame):
        for i in range(1,len(self.region_frames)):
            if self.region_frames[i]==frame:
                self.region_frames[i-1],self.region_frames[i] = self.region_frames[i],self.region_frames[i-1]
                break
        self.update_regionframes_grid()

    def move_frame_down(self,frame):
        for i in range(len(self.region_frames)-1):
            if self.region_frames[i]==frame:
                self.region_frames[i+1],self.region_frames[i] = self.region_frames[i],self.region_frames[i+1]
                break
        self.update_regionframes_grid()

    def upup(self,frame):
        self.region_frames.remove(frame)
        self.region_frames = [frame] + self.region_frames
        self.update_regionframes_grid()
    def downdown(self,frame):
        self.region_frames.remove(frame)
        self.region_frames = self.region_frames + [frame]
        self.update_regionframes_grid()

    def update_regionframes_grid(self):
        for i in range(len(self.region_frames)):
            pady=(2,2)
            if i==0: pady =  (4,pady[1])
            if i==len(self.region_frames)-1: pady =  (pady[0],5)
            self.region_frames[i].grid(row=i+1, column=0, padx=10, pady=pady, sticky="new",columnspan=3)

    def reload(self,params_list=[]):
        n_previous_regions = len(self.region_frames)
        n_new_regions = len(params_list)
        if n_previous_regions>n_new_regions:
            for i in range(n_previous_regions-n_new_regions):
                region_frame = self.region_frames.pop()
                region_frame.grid_forget()
                self.unused_region_frames.append(region_frame)
        for i in range(min(n_new_regions,n_previous_regions)):
            self.region_frames[i].reload(params_list[i])
        if n_new_regions>n_previous_regions:
            for i in range(n_previous_regions,n_new_regions):
                self.add_region_frame(params_list[i],update_grid=False)
        self.update_regionframes_grid()


    def get_regions(self):
        l=[]
        for x in self.region_frames:
           l.append(x.get_region())
        return l

class HighlightFrame(ctk.CTkFrame):
    def __init__(self, master,params={}):
        super().__init__(master,fg_color=colorB,border_width=2)

        with resources.as_file(resources.files(figeno.GUI.data) / "up.png") as infile:
            image_up = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,15))
        with resources.as_file(resources.files(figeno.GUI.data) / "down.png") as infile:
            image_down = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,15))
        self.up_button=ctk.CTkButton(self,width=5,height=5,text="",image=image_up,command=self.up,font=font,fg_color=colorA,hover_color=colorC)
        self.down_button=ctk.CTkButton(self,5,5,text="",image=image_down,command=self.down,font=font,fg_color=colorA,hover_color=colorC)
        self.up_button.grid(row=0,column=0,padx=2,pady=2)
        self.down_button.grid(row=1,column=0,padx=2,pady=2)

        with resources.as_file(resources.files(figeno.GUI.data) / "upup.png") as infile:
            image_upup =ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,15))
        with resources.as_file(resources.files(figeno.GUI.data) / "downdown.png") as infile:
            image_downdown = ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS),size=(15,15))
        self.upup_button=ctk.CTkButton(self,width=5,height=5,text="",image=image_upup,command=self.upup,font=font,fg_color=colorA,hover_color=colorC)
        self.downdown_button=ctk.CTkButton(self,5,5,text="",image=image_downdown,command=self.downdown,font=font,fg_color=colorA,hover_color=colorC)
        self.upup_button.grid(row=0,column=1,padx=2,pady=2)
        self.downdown_button.grid(row=1,column=1,padx=2,pady=2)

        self.chr_label=ctk.CTkLabel(self,10,20,text="chr:",font=font)
        self.chr_label.grid(row=0,column=2,rowspan=2,padx=(10,0))
        self.chr_stringvar = tkinter.StringVar(self, "")
        if "chr" in params: self.chr_stringvar.set(params["chr"])
        self.chr_entry=ctk.CTkEntry(self,50,20,font=font,textvariable=self.chr_stringvar)
        self.chr_entry.grid(row=0,column=3,rowspan=2,padx=(0,10))

        self.start_label=ctk.CTkLabel(self,10,20,text="start:",font=font)
        self.start_label.grid(row=0,column=4,rowspan=2,padx=(10,0))
        self.start_stringvar = tkinter.StringVar(self, "")
        if "start" in params: self.start_stringvar.set(params["start"])
        self.start_entry=ctk.CTkEntry(self,200,20,font=font,textvariable=self.start_stringvar)
        self.start_entry.grid(row=0,column=5,rowspan=2,padx=(0,10))

        self.end_label=ctk.CTkLabel(self,10,10,text="end:",font=font)
        self.end_label.grid(row=0,column=6,rowspan=2,padx=(10,0))
        self.end_stringvar = tkinter.StringVar(self, "")
        if "end" in params: self.end_stringvar.set(str(params["end"]))
        self.end_entry=ctk.CTkEntry(self,200,20,font=font,textvariable=self.end_stringvar)
        self.end_entry.grid(row=0,column=7,rowspan=2,padx=(0,10))

        self.colorbutton=ColorButton(self,25,25,color="#eba434")
        if "color" in params:
            self.colorbutton.set(params["color"])
        self.colorbutton.grid(row=0,column=8,padx=(15,2),rowspan=2)

        self.alpha_label=ctk.CTkLabel(self,10,10,text="Opacity:",font=font)
        self.alpha_label.grid(row=0,column=9,rowspan=2,padx=(10,0))
        self.alpha_stringvar = tkinter.StringVar(self, "0.3")
        if "alpha" in params: self.alpha_stringvar.set(str(params["alpha"]))
        self.alpha_entry=ctk.CTkEntry(self,60,20,font=font,textvariable=self.alpha_stringvar)
        self.alpha_entry.grid(row=0,column=10,rowspan=2,padx=(0,10))
       
        with resources.as_file(resources.files(figeno.GUI.data) / "delete.png") as infile:
            image_delete =ctk.CTkImage(Image.open(infile).resize((20,20),Image.Resampling.LANCZOS))
        self.delete_button=ctk.CTkButton(self,10,10,text="",image=image_delete,command=self.removeframe,font=font,fg_color="#e74c3c",hover_color=color2hovercolor("#e74c3c"))
        self.delete_button.grid(row=0,column=12,rowspan=2,padx=10,sticky="e")
        self.grid_columnconfigure(11,weight=1)

    def removeframe(self):
        self.master.remove_region_frame(self)
    def up(self):
        self.master.move_frame_up(self)
    def down(self):
        self.master.move_frame_down(self)
    def upup(self):
        self.master.upup(self)
    def downdown(self):
        self.master.downdown(self)
    
    def reload(self,params={}):
        if "chr" in params: self.chr_stringvar.set(params["chr"])
        else: self.chr_stringvar.set("")
        if "start" in params: self.start_stringvar.set(str(params["start"]))
        else: self.start_stringvar.set("")
        if "end" in params: self.end_stringvar.set(str(params["end"]))
        else: self.end_stringvar.set("")
        if "color" in params: self.colorbutton.set(params["color"])
        else: self.colorbutton.set("#eba434")
        if "alpha" in params: self.alpha_stringvar.set(str(params["alpha"]))
        else: self.alpha_stringvar.set("0.3")

    def get_region(self):
        params =  {"chr":self.chr_entry.get()}
        try: params["start"]=int(self.start_entry.get())
        except: pass
        try: params["end"]=int(self.end_entry.get())
        except: pass
        params["color"]=self.colorbutton.get()
        params["alpha"] = float(self.alpha_stringvar.get())
        return params
    

