import traceback
import webbrowser
import os
import  sys
import logging
import subprocess
from flask import Flask, jsonify, request
import json
from figeno import figeno_make
from figeno.genes import find_genecoord_wrapper
from figeno.utils import KnownException

app = Flask(__name__, static_folder='./build', static_url_path='/')

if sys.platform=="win32":
    last_dir=os.path.expanduser("~")
else:
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
        #os.startfile(data["file"])
        #webbrowser.open_new_tab("file:///"+data["file"])
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
    debug=True
    if args is not None: debug=args.debug
    if not debug:
        logging.getLogger('werkzeug').disabled = True
    port=5000
    if args is not None: port = args.port
    #print("Starting server on http://localhost:"+str(port))
    #webbrowser.open_new_tab('http://localhost:'+str(port)+'/')
    if args.host!="":
        print("Starting server on "+args.host+":"+str(port))
        app.run(debug=debug,port=port,host=args.host)
    else:
        print("Starting server on http://localhost:"+str(port))
        app.run(debug=debug,port=port)

if __name__=="__main__":
    main()
    
