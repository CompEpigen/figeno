import traceback
import webbrowser
import os
import logging
from flask import Flask, jsonify, request
import json
import crossfiledialog
from figeno import figeno_make

app = Flask(__name__, static_folder='./build', static_url_path='/')




@app.route('/browse')
def browse():
    t=crossfiledialog.open_file()
    return jsonify({"path":t})

@app.route('/open_files')
def open_files():
    t=crossfiledialog.open_multiple()
    return jsonify({"files":t})

@app.route('/save')
def save():
    t=crossfiledialog.save_file()
    return jsonify({"path":t})

@app.route('/save_config', methods = ['POST'])
def save_config():
    if request.is_json:
        data = request.get_json()
        filename=crossfiledialog.save_file()
        if len(filename)>0:
            with open(filename,"w") as fp:
                json.dump(data,fp,indent= "\t")
            return jsonify({"path":filename})
        else: return {}

@app.route('/load_config')
def load_config():
    filename=crossfiledialog.open_file(filter="*.json")
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
    
