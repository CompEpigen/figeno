import traceback
import webbrowser
import os
import sys
import logging
from flask import Flask, jsonify, request
import json
import filedialpy
from figeno import figeno_make
from figeno.genes import find_genecoord_refseq_wrapper

app = Flask(__name__, static_folder='./build', static_url_path='/')

last_dir=os.getcwd()
config_dir=last_dir
config_file="config.json"


#@app.route('/browse/', defaults={'path': ''})
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
        t=filedialpy.openFile(initial_dir=start_dir, title="Select file")
        if len(t)>0: last_dir= os.path.dirname(t)
        return jsonify({"path":t})
    elif data["dialog_type"]=="open_files":
        t=filedialpy.openFiles(initial_dir=last_dir,title="Select files")
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
        t=filedialpy.saveFile(initial_dir=start_dir,initial_file=save_filename,title="Select output path for the figure", filter="*.svg *.pdf *.png *.eps *.ps")
        if len(t)>0:last_dir= os.path.dirname(t)
        return jsonify({"path":t})
    elif data["dialog_type"]=="load_config":
        filename=filedialpy.openFile(initial_dir=config_dir,filter="*.json",title="Select config file to load")
        return jsonify({"path":filename})
    elif data["dialog_type"]=="save_config":
        filename=filedialpy.saveFile(initial_dir=config_dir,initial_file=config_file,title="Select path to save the config file",filter="*.json")
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
    debug=True
    if args is not None: debug=args.debug
    if not debug:
        logging.getLogger('werkzeug').disabled = True
    port=5000
    if args is not None: port = args.port
    if sys.platform=="darwin": web_address="http://127.0.0.1:"+str(port)+"/"
    else:  web_address="http://localhost:"+str(port)+"/"
    print("Starting local server on "+web_address)
    webbrowser.open_new_tab(web_address)
    app.run(debug=debug,port=port)

if __name__=="__main__":
    main()
    
