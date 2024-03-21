import React from 'react';
import {v4 as uuid4} from 'uuid';
import { ColorButton } from './ColorButton';
import "./style.css";
import { DndItem } from './DndItem';
import { DndContainer } from './DndContainer';
import { PathEntry } from './Track';


function Bam({trackID, bam, copy_bam,delete_bam, set_value,className,openColorPanel, setFileDialogData,setFileDialogActive}){
    const fileClass=(bam.file.endsWith(".bam"))?"":"unsetPath";

    function setLabel(e,index){
        const labels = [...bam.labels]
        labels[index] = e.target.value
        set_value("labels",labels)
    }
    function setColor(c,index){
        const colors = [...bam.colors]
        colors[index] = c
        set_value("haplotype_colors",colors)
    }
    function setNgroups(){
        if (bam.split_by_haplotype){
            set_value("labels",[bam.labels[0]]);
            set_value("colors",[bam.colors[0]]);
            set_value("split_by_haplotype",false);
        }
        else{
            set_value("labels",[bam.labels[0],""]);
            set_value("colors",[bam.colors[0],"#e67e22"]);
            set_value("split_by_haplotype",true);
        }
    }
    return(
        <DndItem id={bam.id} copy_item={copy_bam} delete_item={delete_bam} className={className}> 
            <div className="trackOption" style={{display:"flex",alignItems:"center"}}>
                <PathEntry id={"file"+bam.id+trackID} label="File:" value={bam.file} set_value={(val) => set_value("file",val)} className={fileClass}  setFileDialogData={setFileDialogData} setFileDialogActive={setFileDialogActive}/>

                <div className='formItem'>
                    <input id={"base"+bam.id+trackID} style={{width:"2em"}} value={bam.base} onChange={(e) => set_value("base",e.target.value)} />
                </div>

                <div className='formItem'>
                    <input id={"mod"+bam.id+trackID} style={{width:"2em"}} value={bam.mod} onChange={(e) => set_value("mod",e.target.value)} />
                </div>
                <div className='formItem'>
                    <label htmlFor={"min_coverage"+bam.id+trackID}>Min coverage:</label>
                    <input id={"min_coverage"+bam.id+trackID} style={{width:"2em"}} value={bam.min_coverage}  onChange={(e) => set_value("min_coverage",e.target.value)} />
                </div>
                <div className='formItem'>
                    <label htmlFor={"linewidth"+bam.id+trackID}>Linewidth:</label>
                    <input id={"linewidth"+bam.id+trackID} style={{width:"2em"}} value={bam.linewidth}  onChange={(e) => set_value("linewidth",e.target.value)} />
                </div>
                <div className='formItem'>
                    <label htmlFor={"opacity"+bam.id+trackID}>Opacity:</label>
                    <input id={"opacity"+bam.id+trackID} style={{width:"2em"}} value={bam.opacity}  onChange={(e) => set_value("opacity",e.target.value)} />
                </div>
                <div className='formItem'>
                    <label htmlFor={"fix_hardclip"+bam.id+trackID}>Fix hardclip:</label>
                    <input type="checkbox" id={"fix_hardclip"+bam.id+trackID} value={bam.fix_hardclip}  onChange={()=>set_value("fix_hardclip",!bam.fix_hardclip)} />
                </div>
                <div className='formItem'>
                    <label htmlFor={"split_by_haplotype"+bam.id+trackID}>Split by haplotype:</label>
                    <input type="checkbox" id={"split_by_haplotype"+bam.id+trackID} value={bam.split_by_haplotype}  onChange={setNgroups} />
                </div>
                {!bam.split_by_haplotype?(
                    <>
                    <div className='formItem'>
                        <label htmlFor={"color"+bam.id+trackID}>Color:</label>
                        <ColorButton id={"color"+bam.id+trackID} color={bam.colors[0]} setColor={(c)=>setColor(c,0)} openColorPanel={openColorPanel}/>
                    </div>
                    </>
                ):(
                    <>
                    <div className='formItem'>
                        <label htmlFor={"color1"+bam.id+trackID}>Color 1:</label>
                        <ColorButton id={"color1"+bam.id+trackID} color={bam.colors[0]} setColor={(c)=>setColor(c,0)} openColorPanel={openColorPanel}/>
                    </div>
                    <div className='formItem'>
                        <label htmlFor={"color2"+bam.id+trackID}>Color 2:</label>
                        <ColorButton id={"color2"+bam.id+trackID} color={bam.colors[1]} setColor={(c)=>setColor(c,1)} openColorPanel={openColorPanel}/>
                    </div>
                    </>
                )}
                
            </div>
        </DndItem>
    )
}

function Bedmethyl({trackID, bed, copy_bed,delete_bed, set_value,className,openColorPanel, setFileDialogData,setFileDialogActive}){
    const fileClass=(bed.file.endsWith(".bed.gz") || bed.file.endsWith(".bedmethyl.gz") )?"":"unsetPath";
    return(
        <DndItem id={bed.id} copy_item={copy_bed} delete_item={delete_bed} className={className}> 
            <div className="trackOption" style={{display:"flex",alignItems:"center"}}>
                <PathEntry id={"file"+bed.id+trackID} label="File:" value={bed.file} set_value={(val) => set_value("file",val)} className={fileClass}  setFileDialogData={setFileDialogData} setFileDialogActive={setFileDialogActive}/>

                <div className='formItem'>
                    <input id={"mod"+bed.id+trackID} style={{width:"2em"}} value={bed.mod} onChange={(e) => set_value("mod",e.target.value)} />
                </div>
                <div className='formItem'>
                    <label htmlFor={"min_coverage"+bed.id+trackID}>Min coverage:</label>
                    <input id={"min_coverage"+bed.id+trackID} style={{width:"2em"}} value={bed.min_coverage}  onChange={(e) => set_value("min_coverage",e.target.value)} />
                </div>
                <div className='formItem'>
                    <label htmlFor={"linewidth"+bed.id+trackID}>Linewidth:</label>
                    <input id={"linewidth"+bed.id+trackID} style={{width:"2em"}} value={bed.linewidth}  onChange={(e) => set_value("linewidth",e.target.value)} />
                </div>
                <div className='formItem'>
                    <label htmlFor={"opacity"+bed.id+trackID}>Opacity:</label>
                    <input id={"opacity"+bed.id+trackID} style={{width:"2em"}} value={bed.opacity}  onChange={(e) => set_value("opacity",e.target.value)} />
                </div>
                    <div className='formItem'>
                        <label htmlFor={"color"+bed.id+trackID}>Color:</label>
                        <ColorButton id={"color"+bed.id+trackID} color={bed.color} setColor={(c)=>set_value("color",c)} openColorPanel={openColorPanel}/>
                    </div>
                
            </div>
        </DndItem>
    )
}

export function BasemodfreqTrack({track,set_value,openColorPanel, setFileDialogData,setFileDialogActive}){

    function add_bam(){
        set_value("bams",[...track.bams,{id:uuid4(),file:"",base:"C",mod:"m",min_coverage:6,linewidth:3,opacity:1.0,fix_hardclip:false,split_by_haplotype:false,
    labels:[""],colors:["#27ae60"]}])
    }
    function delete_bam(id){
        set_value("bams",track.bams.filter(b=>b.id!=id));
    }
    function copy_bam(id){
        const index = track.bams.findIndex((t)=>t.id==id);
        const newBam = {...track.bams[index],id:uuid4()};
        const newList = [...track.bams]
        newList.splice(index,0,newBam);
        set_value("bams",newList);
    }
    function getBamById (id) {
        return track.bams.find((bam) => bam.id === id);
    }

    function set_value_bam(id,attribute,value){
        set_value("bams",track.bams.map((bam)=>{
            if(bam.id!==id) {return bam }
            else {
                bam[attribute]=value;
                return bam
            }
        }))
    }

    function add_bed(){
        set_value("bedmethyls",[...track.bedmethyls,{id:uuid4(),file:"",mod:"m",min_coverage:6,linewidth:3,opacity:1.0,label:"",color:"#27ae60"}])
    }
    function delete_bed(id){
        set_value("bedmethyls",track.bedmethyls.filter(b=>b.id!=id));
    }
    function copy_bed(id){
        const index = track.bedmehtyls.findIndex((t)=>t.id==id);
        const newBed = {...track.bedmethyls[index],id:uuid4()};
        const newList = [...track.bedmethyls]
        newList.splice(index,0,newBed);
        set_value("bedmethyls",newList);
    }
    function getBedById (id) {
        return track.bedmethyls.find((bed) => bed.id === id);
    }

    function set_value_bed(id,attribute,value){
        set_value("bedmethyls",track.bedmethyls.map((bed)=>{
            if(bed.id!==id) {return bed }
            else {
                bed[attribute]=value;
                return bed
            }
        }))
    }

    const header=<div className="subTrackContainerHeader">
        Bams
    <button className="subButtonHeader" onClick={add_bam}>
    <svg xmlns="http://www.w3.org/2000/svg" height="24" viewBox="0 -960 960 960" width="24"><path d="M440-440H200v-80h240v-240h80v240h240v80H520v240h-80v-240Z" fill="#000000" /></svg>
        Add bam</button>
    </div>;
    function show_active_item(active_id){
        return <Bam bam={{...getBamById(active_id),id:active_id+"overlay"}} className={"track trackOverlay"} /> 
    }

    const header_bed=<div className="subTrackContainerHeader">
    Bedmethyls
    <button className="subButtonHeader" onClick={add_bed}>
    <svg xmlns="http://www.w3.org/2000/svg" height="24" viewBox="0 -960 960 960" width="24"><path d="M440-440H200v-80h240v-240h80v240h240v80H520v240h-80v-240Z" fill="#000000" /></svg>
    Add bedmethyl</button>
    </div>;
    function show_active_item_bed(active_id){
        return <Bedmethyl bed={{...getBedById(active_id),id:active_id+"overlay"}} className={"track trackOverlay"} /> 
    }

    
    return (
        <>
           
            <DndContainer header={header} items={track.bams} setItems={(x)=>set_value("bams",x)} show_active_item={show_active_item} className="subContainer">
            {track.bams.map((bam)=>{
                    return <Bam key={bam.id} trackID={track.id} bam={bam} set_value={(attribute,value)=>set_value_bam(bam.id,attribute,value)}  className={"track noShadow"} copy_bam={()=>copy_bam(bam.id)} delete_bam={()=>delete_bam(bam.id)} openColorPanel={openColorPanel} setFileDialogData={setFileDialogData} setFileDialogActive={setFileDialogActive} />
                })}

            </DndContainer>


           <DndContainer header={header_bed} items={track.bedmethyls} setItems={(x)=>set_value("bedmethyls",x)} show_active_item={show_active_item_bed} className="subContainer subContainerNarrow">
           {track.bedmethyls.map((bed)=>{
                   return <Bedmethyl key={bed.id} trackID={track.id} bed={bed} set_value={(attribute,value)=>set_value_bed(bed.id,attribute,value)}  className={"track noShadow"} copy_bed={()=>copy_bed(bed.id)} delete_bed={()=>delete_bed(bed.id)} openColorPanel={openColorPanel} setFileDialogData={setFileDialogData} setFileDialogActive={setFileDialogActive} />
               })}

           </DndContainer>


        </>
    )
}
