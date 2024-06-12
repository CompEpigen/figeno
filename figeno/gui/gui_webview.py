import traceback
import logging
import os
import sys
import subprocess
from flask import Flask, jsonify, request
import json

from figeno import figeno_make
from figeno.genes import find_genecoord_wrapper
from figeno.utils import KnownException
import webview
if "http_proxy" in os.environ: del os.environ["http_proxy"]
if "HTTP_PROXY" in os.environ: del os.environ["HTTP_PROXY"]
if sys.platform=="win32":
    last_dir=os.path.expanduser("~")
else:
    last_dir=os.getcwd()
config_dir=last_dir
config_file = "config.json"

app = Flask(__name__, static_folder='./build', static_url_path='/')
window = webview.create_window('figeno',app,width=webview.screens[0].width,height=webview.screens[0].height)

@app.route('/browse', methods=["POST"])
def browse():
    global last_dir
    global config_dir
    global config_file
    data=request.get_json()
    if data["dialog_type"]=="open_file":
        start_dir = last_dir
        if len(data["path"])>0 and os.path.exists(os.path.dirname(data["path"])):
            start_dir = os.path.dirname(data["path"])
        t = window.create_file_dialog(webview.OPEN_DIALOG, allow_multiple=False ,directory=start_dir)
        if t is None: t=""
        if not isinstance(t,str):  t=t[0]
        if len(t)>0: last_dir= os.path.dirname(t)
        return jsonify({"path":t})
    elif data["dialog_type"]=="open_files":
        t=window.create_file_dialog(webview.OPEN_DIALOG, allow_multiple=True ,directory=last_dir)
        if t is None: t=[]
        if len(t)>0 and len(t[0])>0: last_dir= os.path.dirname(t[0])
        return jsonify({"files":t})
    elif data["dialog_type"]=="save_file":
        start_dir = last_dir
        if len(data["path"])>0 and os.path.exists(os.path.dirname(data["path"])):
            start_dir=os.path.dirname(data["path"])
        save_filename="figure.svg"
        if len(data["path"])>0:
            filename = os.path.basename(data["path"])
            if len(filename)>0 and (filename.endswith(".svg") or filename.endswith(".pdf") or filename.endswith(".ps") or filename.endswith(".eps") or filename.endswith(".png")):
                save_filename=filename
        t=window.create_file_dialog( webview.SAVE_DIALOG,  save_filename=save_filename,directory=start_dir)
        if t is None: t=""
        if not isinstance(t,str):  t=t[0]
        if len(t)>0:last_dir= os.path.dirname(t)
        return jsonify({"path":t})
    elif data["dialog_type"]=="load_config":
        filename=filename = window.create_file_dialog(webview.OPEN_DIALOG, allow_multiple=False ,directory=config_dir,file_types=('JSON files (*.json)', 'All files (*.*)'))
        if filename is None: filename=""
        if not isinstance(filename,str):  filename=filename[0]
        return jsonify({"path":filename})
    elif data["dialog_type"]=="save_config":
        filename= window.create_file_dialog( webview.SAVE_DIALOG,  save_filename=config_file,directory=config_dir,file_types=('JSON files (*.json)', 'All files (*.*)'))
        if filename is None: filename=""
        if not isinstance(filename,str):  filename=filename[0]
        return jsonify({"path":filename})
@app.route('/save_config', methods = ['POST'])
def save_config():
    global config_dir
    global config_file
    if request.is_json:
        data = request.get_json()
        filename=data["path"]
        if len(filename)>0:
            config_dir=os.path.dirname(filename)
            config_file=os.path.basename(filename)
            with open(filename,"w") as fp:
                json.dump(data["config"],fp,indent= "\t")
            return jsonify({"path":filename})
        else: return {}

@app.route('/load_config', methods=['POST'])
def load_config():
    global config_dir
    global config_file
    data=request.get_json()
    if len(data["path"])>0:
        config_dir=os.path.dirname(data["path"])
        config_file=os.path.basename(data["path"])
        with open(data["path"],"r") as fp:
            config = json.load(fp)
        return jsonify(config)
    else:
        return {}

@app.route('/run', methods = ['POST'])
def run():
    if request.is_json:
        data = request.get_json()
        try:
            warnings=[]
            figeno_make(data,warnings=warnings)
            warning="\n\n".join(warnings)
            #if warning!="": print(warning)
            return jsonify({"status":"success","message":data["output"]["file"],"warning":warning})
        except KnownException as e:
            print(str(e))
            return jsonify({"status":"known_error","message":str(e)})
        except Exception as e:
            print(traceback.format_exc())
            return jsonify({"status":"unknown_error","message":traceback.format_exc()})
        
@app.route('/open_image', methods = ['POST'])
def open_image():
    if request.is_json:
        data = request.get_json()
        if sys.platform == "win32":
            os.startfile(data["file"])
        else:
            opener = "open" if sys.platform == "darwin" else "xdg-open"
            subprocess.call([opener,data["file"]])
    return {}

@app.route('/open_dir', methods = ['POST'])
def open_dir():
    if request.is_json:
        data = request.get_json()
        if sys.platform == "win32":
            os.startfile(os.path.dirname(data["file"]))
        else:
            opener = "open" if sys.platform == "darwin" else "xdg-open"
            subprocess.call([opener,os.path.dirname(data["file"])])
    return {}
        
@app.route('/find_gene', methods = ['POST'])
def find_gene():
    if request.is_json:
        data = request.get_json()
        try:
            chr,start,end=find_genecoord_wrapper(data["gene_name"],data["reference"],data["genes_file"])
            if chr=="": return jsonify({"status":"known_error","message":"Could not find gene: "+data["gene_name"]})
            else: return jsonify({"status":"success","chr":chr,"start":start,"end":end})
        except KnownException as e:
            print(str(e))
            return jsonify({"status":"known_error","message":str(e)})
        except Exception as e:
            print(traceback.format_exc())
            return jsonify({"status":"unknown_error","message":traceback.format_exc()})
        
@app.route('/get_all_chromosomes', methods = ['POST'])
def get_all_chromosomes():
    if request.is_json:
        data = request.get_json()
        try:
            chromosomes=[]
            if "cytobands_file" in data and data["cytobands_file"]!="":
                if not os.path.isfile(data["cytobands_file"]): raise KnownException("The cytobands file could not be found: "+str(data["cytobands_file"]+". Will use the human chromosomes by default."))
                with open(data["cytobands_file"],"r") as infile:
                    for line in infile:
                        if line.startswith("#"): continue
                        linesplit = line.rstrip("\n").split("\t")
                        chr = linesplit[0].lstrip("chr")
                        if not chr in chromosomes: chromosomes.append(chr)
                return jsonify({"status":"success","chromosomes":chromosomes})
            else:
                return jsonify({"status":"known_error","message":"No cytobands (or .fai) file was provided, so by default all human chromosomes were added. Please provide a cytobands (or .fai) file if you want to add the chromosomes present in your reference.","chromosomes":[]})
        except KnownException as e:
            print(str(e))
            return jsonify({"status":"known_error","message":str(e),"chromosomes":[]})
        except Exception as e:
            print(traceback.format_exc())
            return jsonify({"status":"unknown_error","message":traceback.format_exc(),"chromosomes":[]})



@app.route('/')
def ind():
    return app.send_static_file('index.html')

def main(args=None):
    
    debug=False
    if args is not None: debug=args.debug
    if not debug:
        logging.getLogger('werkzeug').disabled = True
    webview.start(debug=args.debug)