import traceback
import logging
import os
from flask import Flask, jsonify, request
import json

from figeno import figeno_make
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
    data=request.get_json()
    start_dir = last_dir
    if len(data["path"])>0 and os.path.exists(os.path.dirname(data["path"])):
        start_dir = os.path.dirname(data["path"])
    result = window.create_file_dialog(
        webview.OPEN_DIALOG, allow_multiple=False ,directory=start_dir)
    
    if result is not None:
        filename = result[0]
        last_dir=os.path.dirname(filename)
        return jsonify({"path":filename})
    else:
        return jsonify({"path":""})

@app.route('/open_files')
def open_files():
    global last_dir
    result = window.create_file_dialog(
        webview.OPEN_DIALOG, allow_multiple=True ,directory=last_dir)
    if result is not None:
        last_dir = os.path.dirname(result[0])

    return jsonify({"files":result})

@app.route('/save',methods=["POST"])
def save():
    global last_dir
    data=request.get_json()
    start_dir = last_dir
    save_filename="figure.svg"
    if len(data["path"])>0 and os.path.exists(os.path.dirname(data["path"])):
        start_dir=os.path.dirname(data["path"])
    if len(data["path"])>0:
        filename = os.path.basename(data["path"])
        if len(filename)>0 and (filename.endswith(".svg") or filename.endswith(".pdf") or filename.endswith(".ps") or filename.endswith(".eps") or filename.endswith(".png")):
            save_filename=filename
    t = window.create_file_dialog( webview.SAVE_DIALOG,  save_filename=save_filename,directory=start_dir)
    if t is not None:
        if not isinstance(t,str):  t=t[0]
        last_dir = os.path.dirname(t)
    return jsonify({"path":t})

@app.route('/save_config', methods = ['POST'])
def save_config():
    global config_dir
    global config_file
    if request.is_json:
        data = request.get_json()
        filename = window.create_file_dialog( webview.SAVE_DIALOG,  save_filename=config_file,directory=config_dir,file_types=('JSON files (*.json)', 'All files (*.*)'))
        if filename is not None:
            if not isinstance(filename,str): filename = filename[0]
            config_dir=os.path.dirname(filename)
            config_file = os.path.basename(filename)
            with open(filename,"w") as fp:
                json.dump(data,fp,indent= "\t")
            return jsonify({"path":filename})
        else:
            return jsonify({})

@app.route('/load_config')
def load_config():
    global config_dir
    global config_file
    
    filename = window.create_file_dialog(
        webview.OPEN_DIALOG, allow_multiple=False ,directory=config_dir,file_types=('JSON files (*.json)', 'All files (*.*)'))
    
    if filename is not None:
        filename = filename[0]
        config_dir=os.path.dirname(filename)
        config_file=os.path.basename(filename)
        with open(filename,"r") as fp:
            config = json.load(fp)
        return jsonify(config)
    else:
        return jsonify({})

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


@app.route('/')
def ind():
    return app.send_static_file('index.html')

def main(args=None):
    
    debug=False
    if args is not None: debug=args.debug
    if not debug:
        logging.getLogger('werkzeug').disabled = True
    webview.start(debug=args.debug)