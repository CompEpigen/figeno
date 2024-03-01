import {React, Fragment,useState,useRef, useCallback} from 'react';
import { useSortable } from '@dnd-kit/sortable';
import {CSS} from '@dnd-kit/utilities';
import { HexColorPicker, HexColorInput } from "react-colorful";
import useClickOutside from "./useClickOutside";
import "./style.css";

export function ColorButton({color,setColor,openColorPanel,id}){
    // selectColor opens the ColorPanel

    //
    return (
    <>
        <button  id={id} className="ColorButton" style={{ backgroundColor: color }} onClick={() => openColorPanel(color,setColor)}/>
    </>
    )
}

export function ColorPanel({color,setColor,updateColor,close}){
    const predefined_colors=["#000000","#222222","#444444","#666666","#888888","#aaaaaa",
                            "#c0392b","#e74c3c","#d35400","#e67e22","#f39c12","#f1c40f",
                            "#2980b9","#3498db","#27ae60","#2ecc71","#16a085","#1abc9c",
                            "#8e44ad","#9b59b6","#2c3e50","#34495e","#a29bfe","#ffffff"]
    //const [color,setColor] = useState("#333333");
    // setColor sets the color of the colorpanel, updateColor updates the color of the color button.
    const popover = useRef();
    useClickOutside(popover, close);

    function handleOK(){
        updateColor(color);
        close();
    }
    return(
    <div className="ColorPanel" style={{backgroundColor:color}} ref={popover}>
        <div style={{backgroundColor:"white",padding:"10px",width:"inherit", borderRadius:"10px", 
        display:"flex",justifyContent:"center",alignContent:"center",alignItems:"center",flexDirection:"column"}}>
                <h2>Predefined colors</h2>
                <div className="PredefinedColorsBlock">
                    {predefined_colors.map((color)=>(
                        <ColorButtonInPanel key={"colorbutton"+color} color={color} setColor={setColor}/>
                    ))}
                </div>
                <h2>Custom colors</h2>
                <div style={{display:"flex",flexDirection:"column", width:"inherit",justifyContent:"center",alignContent:"center",alignItems:"center",gap:"10px"}}>
                    <HexColorPicker color={color} onChange={setColor} />
                    <input style={{fontSize:"1.3em", width:"5em"}} id={"colorpanelInput"} value={color} onChange={(e) => setColor(e.target.value)} ></input>
                    <div style={{display:"flex", gap:"5px"}}>
                        <button style={{fontSize:"1.3em"}} id={"colorpanelCancel"} onClick={close}>Cancel</button>
                        <button style={{fontSize:"1.3em"}}  id={"colorpanelOK"} onClick={handleOK}>OK</button>
                    </div>
                </div>
        </div>
    </ div>
    )
}

function ColorButtonInPanel({color,setColor}){
    return (
    <>
        <button  id={"color"+color} className="ColorButton" style={{ backgroundColor: color }} onClick={() => setColor(color)}/>
    </>
    )
}