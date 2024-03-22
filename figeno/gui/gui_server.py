import traceback
import webbrowser
import os
import logging
from flask import Flask, jsonify, request
import json
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
    data=request.get_json()
    start_dir = last_dir
    if data["dialog_type"] in ["save_config","load_config"]: start_dir=config_dir
    if len(data["path"])>0 and os.path.exists((data["path"])):
        if os.path.isdir(data["path"]):
            start_dir = data["path"]
        else: start_dir = os.path.dirname(data["path"])
        last_dir=start_dir
    dirs=[]
    files=[]
    for f in os.listdir(start_dir):
        full_path = os.path.join(start_dir,f)
        if os.path.isfile(full_path): 
            if (not data["dialog_type"] in ["load_config","save_config"]) or f.endswith(".json"):
                files.append(f)
        else: dirs.append(f)
    return jsonify({"current_dir":start_dir,"dirs":sorted(dirs),"files":sorted(files)})

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
    print("Starting server on http://localhost:"+str(port))
    #webbrowser.open_new_tab('http://localhost:'+str(port)+'/')
    app.run(debug=debug,port=port)

if __name__=="__main__":
    main()
    
