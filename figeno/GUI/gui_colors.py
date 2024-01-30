import tkinter
from tkinter.colorchooser import askcolor
import customtkinter as ctk
import matplotlib.colors as mcolors

font=("arial",16)
colorD=("gray70", "gray25")
colorC=("gray85", "gray16")
colorB=("gray90", "gray13")
colorA=("gray95", "gray10")

def color2rgb(c):
    if c in mcolors.get_named_colors_mapping():
        c = mcolors.get_named_colors_mapping()[c]
    # Assume c is now a hex value
    return mcolors.hex2color(c)

def color2hovercolor(c):
    c = color2rgb(c)
    if c[0]+c[1]+c[2]<0.3: change=+0.3
    else: change=-0.3
    darkened = [max(0,a+change) for a in c]
    return mcolors.rgb2hex(darkened)

class ColorSelector(ctk.CTkToplevel):
    def __init__(self,master,colorbutton,default_color=""):
        super().__init__(master)
        self.attributes('-topmost', True)
        x = master.winfo_x()
        y = master.winfo_y()
        self.geometry("+%d+%d" % (x + 200, y + 200))
        self.minsize(300,300)
        #self.after(1000,self.attributes,'-topmost', False)
        self.title("Color selector")
        self.update()
        self.master = master
        self.colorbutton = colorbutton
        self.color = default_color

        # Predefined frame
        self.predefined_frame = ctk.CTkFrame(self,fg_color=colorC)
        ctk.CTkLabel(self.predefined_frame,text="Predefined colors",font=font).grid(row=0,column=0,columnspan=100,padx=5,pady=5)
        #colors = ["black","#27ae60","#2980b9","#e74c3c","#f1c40f","#2c3e50","#8e44ad","#d35400","#7f8c8d"]
        colors = ["#000000","#222222","#444444","#666666","#888888","#aaaaaa"]
        colors = colors+ ["#c0392b","#e74c3c","#d35400","#e67e22","#f39c12","#f1c40f"]
        colors = colors+ ["#2980b9","#3498db","#27ae60","#2ecc71","#16a085","#1abc9c"]
        colors = colors + ["#8e44ad","#9b59b6","#2c3e50","#34495e","#a29bfe","#ffffff"]
        for i in range(4):
            for j in range(6):
                c = colors[i*6+j]
                button = ctk.CTkButton(self.predefined_frame,width=20,height=20,corner_radius=0,fg_color=c,command = self.create_color_func(c),hover_color=color2hovercolor(c),text="") # command=self.create_color_func(colors[i*3+j])
                button.grid(row=i+1,column=j,padx=1,pady=1)

        # Manual frame
        self.manual_frame=ctk.CTkFrame(self,fg_color=colorC)
        self.label_manual= ctk.CTkLabel(self.manual_frame,text="Manual entry",font=font)
        self.entryvar = tkinter.StringVar(self.manual_frame, self.color)
        self.entry=ctk.CTkEntry(self.manual_frame,textvariable=self.entryvar,width=80,font=font,height=12)
        self.validate_button = ctk.CTkButton(self.manual_frame,20,20,text="OK",font=font,command=self.select_color_manual)
        self.label_manual.grid(row=0,column=0,padx=5,pady=5,columnspan=2)
        self.entry.grid(row=1,column=0,padx=5,pady=5)
        self.validate_button.grid(row=1,column=1,padx=5,pady=5)
        self.manual_frame.grid_columnconfigure((0,1),weight=1)

        # Other frame
        self.other_frame=ctk.CTkFrame(self,fg_color=colorC)
        self.other_frame.grid_columnconfigure(0, weight=1)
        self.other_button = ctk.CTkButton(self.other_frame,text="Select other color",font=font,command=self.color_chooser)
        self.other_button.grid(row=0,column=0,padx=5,pady=5,columnspan=100)

        self.predefined_frame.grid(row=0,column=0,padx=5,pady=5,sticky="")
        self.manual_frame.grid(row=1,column=0,padx=5,pady=5,sticky="we")
        self.other_frame.grid(row=2,column=0,padx=5,pady=5,sticky="we")
        self.grid_columnconfigure(0,weight=1)

    def create_color_func(self,color):
        def color_func():
            self.color = color
            self.destroy()
        return color_func
    def select_color_manual(self):
        self.color = self.entryvar.get()
        self.destroy()
    def color_chooser(self):
        rgb,hex = askcolor(parent=self)
        if hex is not None:
            self.color = hex
            self.destroy()
    def destroy(self):
        self.colorbutton.configure(fg_color=self.color)
        self.colorbutton.configure(hover_color=color2hovercolor(self.color))
        super().destroy()

class ColorButton(ctk.CTkButton):
    def __init__(self,master,width=20,height=20,color="#555555",hover_color=None,border_width=1):
        if hover_color is None: hover_color = color2hovercolor(color)
        super().__init__(master,width=width,height=height,fg_color=color,hover_color=hover_color,border_width=border_width,text="",command=self.select_color)

    def select_color(self):
        colorselector = ColorSelector(self,self,default_color=self._fg_color)

    def set(self,color,hover_color=None):
        if hover_color is None: hover_color=color2hovercolor(color)
        self.configure(fg_color=color)
        self.configure(hover_color=hover_color)
    
    def get(self):
        return self.cget("fg_color")
