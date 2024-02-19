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
    //const [color,setColor] = useState("#333333");
    // setColor sets the color of the colorpanel, updateColor updates the color of the color button.
    const popover = useRef();
    useClickOutside(popover, close);

    function handleOK(){
        updateColor(color);
        close();
    }
    return(
    <div className="ColorPanel" ref={popover}>
                <HexColorPicker color={color} onChange={setColor} />
                <input id={"colorpanelInput"} value={color} onChange={(e) => setColor(e.target.value)} ></input>
                <button id={"colorpanelCancel"} onClick={close}>Cancel</button>
                <button id={"colorpanelOK"} onClick={handleOK}>OK</button>
    </ div>
    )
}
