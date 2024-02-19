import traceback
import webbrowser
import os
import logging
from flask import Flask, jsonify, request
import json
import tkinter.filedialog
from tkinter import Tk
from figeno import figeno_make

app = Flask(__name__, static_folder='./build', static_url_path='/')




@app.route('/browse')
def browse():
    root = Tk()
    root.withdraw()
    t=tkinter.filedialog.askopenfilename()
    root.destroy()
    return jsonify({"path":t})

@app.route('/open_files')
def open_files():
    root = Tk()
    root.withdraw()
    t=tkinter.filedialog.askopenfilenames()
    root.destroy()
    return jsonify({"files":t})

@app.route('/save')
def save():
    root = Tk()
    root.withdraw()
    t=tkinter.filedialog.asksaveasfilename()
    root.destroy()
    return jsonify({"path":t})

@app.route('/save_config', methods = ['POST'])
def save_config():
    if request.is_json:
        data = request.get_json()
        root = Tk()
        root.withdraw()
        filename=tkinter.filedialog.asksaveasfilename(initialfile="config.json")
        #print(t)
        root.destroy()
        print(data)
        with open(filename,"w") as fp:
            json.dump(data,fp,indent= "\t")
        return jsonify({"path":filename})

@app.route('/load_config')
def load_config():
    root = Tk()
    root.withdraw()
    filename=tkinter.filedialog.askopenfilename()
    root.destroy()
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
    debug=True
    if args is not None: debug=args.debug
    if not debug:
        logging.getLogger('werkzeug').disabled = True
    port=5000
    if args is not None: port = args.port
    
    webbrowser.open_new_tab('http://localhost:'+str(port)+'/')
    app.run(debug=debug,port=port)

if __name__=="__main__":
    main()
    
