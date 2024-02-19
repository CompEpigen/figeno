import React from 'react';
import { ColorButton } from './ColorButton';
import "./style.css";
import { DndItem } from './DndItem';





export function Highlight({region, set_value, className, copy_region,delete_region, openColorPanel}) {
    const chrClass=(region.chr!=="")?"":"unsetPath";
    return (
        <DndItem id ={region.id} copy_item={()=>copy_region(region.id)} delete_item={()=>delete_region(region.id)} className={className}>
            <div className="trackHeader">
                    
                <div className='formItem'>
                    <label htmlFor={"chr"+region.id}>chr:</label>
                    <input id={"chr"+region.id} value ={region.chr} style={{width:"2em"}} onChange={(e) => set_value("chr",e.target.value)} className={chrClass} />
                </div>

                <div className='formItem'>
                    <label htmlFor={"start"+region.id}>start:</label>
                    <input id={"start"+region.id} value ={region.start} style={{width:"6em"}} onChange={(e) => set_value("start",e.target.value)}></input>
                </div>

                <div className='formItem'>
                    <label htmlFor={"end"+region.id}>end:</label>
                    <input id={"end"+region.id} value ={region.end} style={{width:"6em"}} onChange={(e) => set_value("end",e.target.value)}></input>
                </div>
                <ColorButton id={"color"+region.id} color={region.color} setColor={(c)=>set_value("color",c)} openColorPanel={openColorPanel}/>
                
                <div className='formItem'>
                    <label htmlFor={"opacity"+region.id}>Opacity:</label>
                    <input id={"opacity"+region.id} value ={region.opacity} style={{width:"2em"}} onChange={(e) => set_value("opacity",e.target.value)}></input>
                </div>
            </div>
        </DndItem>
    )
}
