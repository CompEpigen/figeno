import traceback
import logging
from flask import Flask, jsonify, request
import json

from figeno import figeno_make
import webview

app = Flask(__name__, static_folder='./build', static_url_path='/')
window = webview.create_window('figeno',app,width=webview.screens[0].width,height=webview.screens[0].height)



@app.route('/browse')
def browse():
    result = window.create_file_dialog(
        webview.OPEN_DIALOG, allow_multiple=False )
    if result is not None:
        return jsonify({"path":result[0]})
    else:
        return jsonify({"path":""})

@app.route('/open_files')
def open_files():
    result = window.create_file_dialog(
        webview.OPEN_DIALOG, allow_multiple=True )

    return jsonify({"files":result})

@app.route('/save')
def save():
    t = window.create_file_dialog( webview.SAVE_DIALOG, directory='/', save_filename='figure.svg')
    return jsonify({"path":t})

@app.route('/save_config', methods = ['POST'])
def save_config():
    if request.is_json:
        data = request.get_json()
        filename = window.create_file_dialog( webview.SAVE_DIALOG, directory='/', save_filename='config.json')
        #filename=tkinter.filedialog.asksaveasfilename(initialfile="config.json")
        #print(t)
        with open(filename,"w") as fp:
            json.dump(data,fp,indent= "\t")
        return jsonify({"path":filename})

@app.route('/load_config')
def load_config():
    
    filename = window.create_file_dialog(
        webview.OPEN_DIALOG, allow_multiple=False )[0]
    
    if len(filename)>0:
        with open(filename,"r") as fp:
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


@app.route('/')
def ind():
    return app.send_static_file('index.html')

def main(args=None):
    
    debug=False
    if args is not None: debug=args.debug
    if not debug:
        logging.getLogger('werkzeug').disabled = True
    webview.start(debug=args.debug)