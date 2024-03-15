import traceback
import webbrowser
import os
import logging
from flask import Flask, jsonify, request
import json
import filedialpy
from figeno import figeno_make

app = Flask(__name__, static_folder='./build', static_url_path='/')

last_dir=os.getcwd()
config_dir=last_dir
config_file="config.json"


#@app.route('/browse/', defaults={'path': ''})
@app.route('/browse', methods=["POST"])
def browse():
    global last_dir
    data=request.get_json()
    start_dir = last_dir
    if len(data["path"])>0 and os.path.exists(os.path.dirname(data["path"])):
        start_dir = os.path.dirname(data["path"], title="Select file")
    t=filedialpy.openFile(initial_dir=start_dir)
    if len(t)>0: last_dir= os.path.dirname(t)
    return jsonify({"path":t})

@app.route('/open_files')
def open_files():
    global last_dir
    t=filedialpy.openFiles(initial_dir=last_dir,title="Select files")
    if len(t)>0 and len(t[0])>0: last_dir= os.path.dirname(t[0])
    return jsonify({"files":t})

@app.route('/save',methods=["POST"])
def save():
    global last_dir
    data=request.get_json()
    start_dir = last_dir
    if len(data["path"])>0 and os.path.exists(os.path.dirname(data["path"])):
        start_dir=os.path.dirname(data["path"])
    save_filename="figure.svg"
    if len(data["path"])>0:
        filename = os.path.basename(data["path"])
        if len(filename)>0 and (filename.endswith(".svg") or filename.endswith(".pdf") or filename.endswith(".ps") or filename.endswith(".eps") or filename.endswith(".png")):
            save_filename=filename
    t=filedialpy.saveFile(initial_dir=start_dir,initial_file=save_filename,title="Select output path for the figure", filter="*.svg *.pdf *.png *.eps *.ps")
    if len(t)>0:last_dir= os.path.dirname(t)
    return jsonify({"path":t})

@app.route('/save_config', methods = ['POST'])
def save_config():
    global config_dir
    global config_file
    if request.is_json:
        data = request.get_json()
        filename=filedialpy.saveFile(initial_dir=config_dir,initial_file=config_file,title="Select path to save the config file",filter="*.json")
        if len(filename)>0:
            config_dir=os.path.dirname(filename)
            config_file=os.path.basename(filename)
            with open(filename,"w") as fp:
                json.dump(data,fp,indent= "\t")
            return jsonify({"path":filename})
        else: return {}

@app.route('/load_config')
def load_config():
    global config_dir
    global config_file
    filename=filedialpy.openFile(initial_dir=config_dir,filter="*.json",title="Select config file to load")
    if len(filename)>0:
        config_dir=os.path.dirname(filename)
        config_file=os.path.basename(filename)
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
    print("Starting local server on http://localhost:"+str(port))
    webbrowser.open_new_tab('http://localhost:'+str(port)+'/')
    app.run(debug=debug,port=port)

if __name__=="__main__":
    main()
    
