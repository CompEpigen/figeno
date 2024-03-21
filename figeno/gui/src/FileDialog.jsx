import {React, Fragment} from 'react';
import "./style.css";



export function FileDialog({fileDialogData, setFileDialogData,close}){

    function update_data(new_dir){
        fetch('/browse',{headers: {'Content-Type': 'application/json'}, body: JSON.stringify({path:new_dir,dialog_type:fileDialogData.dialog_type}),
        method:"POST"}).then(res => res.json()).then(data => {
            setFileDialogData({...fileDialogData,current_dir:data.current_dir,dirs:data.dirs,files:data.files})
        });
    }
    function handle_click(dirname){
        let new_dir=fileDialogData.current_dir+"/"+dirname;
        if (fileDialogData.current_dir.slice(-1)=="/") new_dir=fileDialogData.current_dir+dirname;
        update_data(new_dir); 
    }
    function handle_click_file(filename){
        let path=fileDialogData.current_dir+"/"+filename;
        if (fileDialogData.current_dir.slice(-1)=="/") path=fileDialogData.current_dir+filename;
        if (fileDialogData.dialog_type=="save_file" || fileDialogData.dialog_type=="save_config"){
            setFileDialogData({...fileDialogData,file:filename})
        }
        else{
            fileDialogData.update(path);
            close();
        }
        
    }

    function handle_up(){
        let s=fileDialogData.current_dir;
        if (s.slice(-1)=="/") s=s.substring(0,s.length-1);
        s=s.substring(0,s.lastIndexOf("/"));
        update_data(s);
    }
    function handle_refresh(){
        update_data(fileDialogData.current_dir);
    }

    function handle_save(){
        let path=fileDialogData.current_dir+"/"+fileDialogData.file;
        if (fileDialogData.current_dir.slice(-1)=="/") path=fileDialogData.current_dir+fileDialogData.file;
        fileDialogData.update(path);
        close()
    }

    return(
    <div className="LoadingScreen">
        <div style={{gap:"10px",display:"flex"}}>
            <button onClick={handle_up}>
            <svg xmlns="http://www.w3.org/2000/svg" height="24" viewBox="0 -960 960 960" width="24"><path d="M440-160v-487L216-423l-56-57 320-320 320 320-56 57-224-224v487h-80Z"/></svg>
            </button>
            <input style={{fontSize:"1.3em", width:"25em"}} id={"FileDialogDir"} value={fileDialogData.current_dir} onChange={(e) => setFileDialogData({...fileDialogData,current_dir:e.target.value})} />
            <button onClick={handle_refresh}>
            <svg xmlns="http://www.w3.org/2000/svg" height="24" viewBox="0 -960 960 960" width="24"><path d="M480-160q-134 0-227-93t-93-227q0-134 93-227t227-93q69 0 132 28.5T720-690v-110h80v280H520v-80h168q-32-56-87.5-88T480-720q-100 0-170 70t-70 170q0 100 70 170t170 70q77 0 139-44t87-116h84q-28 106-114 173t-196 67Z"/></svg>
            </button>
        </div>

        <div style={{overflowY: "scroll", maxHeight:"600px", display:"flex",flexDirection:"column",margin:"5px"}}>
            {fileDialogData.dirs.map((directory)=>{return <DirectoryButton key={"fileDialog"+directory} dirname={directory} onClick={()=>handle_click(directory)}/>})}
            {fileDialogData.files.map((file)=>{return <FileButton key={"fileDialog"+file} filename={file} onClick={()=>handle_click_file(file)}/>})}
        </div>
        {(fileDialogData.dialog_type=="save_config" || fileDialogData.dialog_type=="save_file")?(
            <div>
                <input style={{fontSize:"1.0em", width:"10em", margin:"5px"}} id={"saveFilename"} value={fileDialogData.file} onChange={(e) => setFileDialogData({...fileDialogData,file:e.target.value})} />
                <button onClick={handle_save}>OK</button>
            </div>
        ):""}
        <button onClick={close}>Close</button>
    </ div>
    )
}

function DirectoryButton({dirname,onClick}){
    return(
        <button className="file" onClick={onClick}>
            <svg xmlns="http://www.w3.org/2000/svg" height="24" viewBox="0 -960 960 960" width="24"><path d="M160-160q-33 0-56.5-23.5T80-240v-480q0-33 23.5-56.5T160-800h240l80 80h320q33 0 56.5 23.5T880-640v400q0 33-23.5 56.5T800-160H160Zm0-80h640v-400H447l-80-80H160v480Zm0 0v-480 480Z"/></svg>
            <div style={{margin:"5px"}}>{dirname}</div>
        </button>
    )
}

function FileButton({filename,onClick}){
    return(
        <button className="file" onClick={onClick}>
            <div style={{margin:"5px"}}>{filename}</div>
        </button>
    )
}
