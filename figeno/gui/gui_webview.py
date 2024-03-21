import traceback
import logging
import os
from flask import Flask, jsonify, request
import json

from figeno import figeno_make
from figeno.genes import find_genecoord_refseq_wrapper
import webview

if "http_proxy" in os.environ: del os.environ["http_proxy"]
if "HTTP_PROXY" in os.environ: del os.environ["HTTP_PROXY"]
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
    print(data)
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
            figeno_make(data)
            return jsonify({"success":True})
        except Exception as e:
            print(traceback.format_exc())
            return jsonify({"success":False,"error":traceback.format_exc()})
        
@app.route('/find_gene', methods = ['POST'])
def find_gene():
    if request.is_json:
        data = request.get_json()
        try:
            chr,start,end=find_genecoord_refseq_wrapper(data["gene_name"],data["reference"],data["genes_file"])
            if chr=="": return jsonify({"success":False,"error":"Could not find gene: "+data["gene_name"]})
            else: return jsonify({"success":True,"chr":chr,"start":start,"end":end})
        except Exception as e:
            print(traceback.format_exc())
            return jsonify({"success":False,"error":traceback.format_exc()})


@app.route('/')
def ind():
    return app.send_static_file('index.html')

def main(args=None):
    
    debug=False
    if args is not None: debug=args.debug
    if not debug:
        logging.getLogger('werkzeug').disabled = True
    webview.start(debug=args.debug)