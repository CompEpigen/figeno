import React from 'react';
import "./style.css";


export function OutputContainer({outputParams,setOutputParams, setFileDialogData, setFileDialogActive}) {

    const fileClass=(outputParams.file.endsWith(".svg") || outputParams.file.endsWith(".pdf") || outputParams.file.endsWith(".eps") 
                        || outputParams.file.endsWith(".ps") || outputParams.file.endsWith(".png"))? "":"unsetPath";

    function set_value(attribute,value){
        const newVal={...outputParams};
        newVal[attribute]= value;
        setOutputParams(newVal);
    }
    

    const header=<div className="TrackContainerHeader">
        <h3>Output</h3>
    </div>;

    return (
        <div className="TrackContainer">
            {header}
            <div className="TrackContainerContent">
                <div className="trackOption" style={{paddingLeft:"15px"}}>
                    <PathEntrySave id={"file"} label="File:" value={outputParams.file} set_value={(val) => set_value("file",val)}  setFileDialogData={setFileDialogData} setFileDialogActive={setFileDialogActive} className={fileClass}/>

                    <div className='formItem'>
                        <label htmlFor={"dpi"}>dpi:</label>
                        <input id={"dpi"} style={{width:"2em"}}  value={outputParams.dpi}  onChange={(e) => set_value("dpi",e.target.value)} ></input>
                    </div>

                    <div className='formItem'>
                        <label htmlFor={"width"}>Width (mm):</label>
                        <input id={"width"} style={{width:"2em"}}  value={outputParams.width}  onChange={(e) => set_value("width",e.target.value)} ></input>
                    </div>

                </div>
            </div>
        </div>
    )
}

function select_file_save(set_value,value, setFileDialogData, setFileDialogActive){
    fetch('/browse',{headers: {'Content-Type': 'application/json'}, body: JSON.stringify({path:value,dialog_type:"save_file"}),
    method:"POST"}).then(res => res.json()).then(data => {
        if (data.hasOwnProperty("current_dir")){
            setFileDialogData({current_dir:data.current_dir,dirs:data.dirs,files:data.files,file:"figure.svg",dialog_type:"save_file",update:function(f){set_value(f)}})
            setFileDialogActive(true);
        }
        else{
            if (data.path.length>0) set_value(data.path);
        }
      });
}
export function PathEntrySave({id,label,value,set_value,className="" ,setFileDialogData, setFileDialogActive}){
    return(
        <div className='formItem'>
            <label htmlFor={id}>{label}</label>
            <input id={id} value={value}  onChange={(e)=>set_value(e.target.value) } style={{marginLeft:"2mm", marginRight:"1mm"}} className={className}></input>
            <button name="file" className="buttonImage" onClick={()=>select_file_save(set_value,value, setFileDialogData, setFileDialogActive)}>
            <svg width="4.6mm" height="4.5mm" version="1.1" viewBox="0 0 576 512" xmlns="http://www.w3.org/2000/svg">
                <path d="M384 480h48c11.4 0 21.9-6 27.6-15.9l112-192c5.8-9.9 5.8-22.1 .1-32.1S555.5 224 544 224H144c-11.4 0-21.9 6-27.6 15.9L48 357.1V96c0-8.8 7.2-16 16-16H181.5c4.2 0 8.3 1.7 11.3 4.7l26.5 26.5c21 21 49.5 32.8 79.2 32.8H416c8.8 0 16 7.2 16 16v32h48V160c0-35.3-28.7-64-64-64H298.5c-17 0-33.3-6.7-45.3-18.7L226.7 50.7c-12-12-28.3-18.7-45.3-18.7H64C28.7 32 0 60.7 0 96V416c0 35.3 28.7 64 64 64H87.7 384z" fill="#555555"/>
            </svg>
            </button>
        </div>
    )
}