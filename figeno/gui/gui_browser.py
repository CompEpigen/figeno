import traceback
import webbrowser
import os
import logging
from flask import Flask, jsonify, request
import json
import crossfiledialog
from figeno import figeno_make

app = Flask(__name__, static_folder='./build', static_url_path='/')

last_dir=os.path.expanduser("~")
config_dir=last_dir



#@app.route('/browse/', defaults={'path': ''})
@app.route('/browse', methods=["POST"])
def browse():
    global last_dir
    data=request.get_json()
    start_dir = last_dir
    if len(data["path"])>0 and os.path.exists(os.path.dirname(data["path"])):
        start_dir = os.path.dirname(data["path"])
    t=crossfiledialog.open_file(start_dir=start_dir)
    if len(t)>0: last_dir= os.path.dirname(t)
    return jsonify({"path":t})

@app.route('/open_files')
def open_files():
    global last_dir
    t=crossfiledialog.open_multiple(start_dir=last_dir)
    if len(t)>0 and len(t[0])>0: last_dir= os.path.dirname(t[0])
    return jsonify({"files":t})

@app.route('/save',methods=["POST"])
def save():
    global last_dir
    data=request.get_json()
    print(data)
    start_dir = last_dir
    if len(data["path"])>0 and os.path.exists(os.path.dirname(data["path"])):
        start_dir=os.path.dirname(data["path"])
    t=crossfiledialog.save_file(start_dir=start_dir)
    if len(t)>0:last_dir= os.path.dirname(t)
    return jsonify({"path":t})

@app.route('/save_config', methods = ['POST'])
def save_config():
    global config_dir
    if request.is_json:
        data = request.get_json()
        filename=crossfiledialog.save_file(start_dir=config_dir)
        if len(filename)>0:
            config_dir=os.path.dirname(filename)
            with open(filename,"w") as fp:
                json.dump(data,fp,indent= "\t")
            return jsonify({"path":filename})
        else: return {}

@app.route('/load_config')
def load_config():
    global config_dir
    filename=crossfiledialog.open_file(filter="*.json",start_dir=config_dir)
    if len(filename)>0:
        config_dir=os.path.dirname(filename)
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
    print("starting local server on http://localhost:"+str(port))
    webbrowser.open_new_tab('http://localhost:'+str(port)+'/')
    app.run(debug=debug,port=port)

if __name__=="__main__":
    main()
    
