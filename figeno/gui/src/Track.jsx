import React from 'react';
import { useState } from 'react';
import { ColorButton } from './ColorButton';
import "./style.css";
import { DndItem } from './DndItem';
import { DndContainer } from './DndContainer';
import {BasemodfreqTrack} from './Basemod.jsx'


const track_types=["","chr_axis", "genes","bed","bigwig","coverage","alignments","hic","basemod_freq","sv","copynumber"];



export function Track({track, set_value, change_track_type, copy_track, delete_track, className, openColorPanel, openTrackTypePanel}) {

    function handleTypeChange(e){
        if (e.target.value!=track.type){
            set_value("type",e.target.value)
        }
    }
    function optionTrack(){
        if (track.type=="chr_axis"){
            return ChrTrack({"track":track,"set_value":set_value});
        }
        if (track.type=="genes"){
            return GenesTrack({"track":track,"set_value":set_value,"openColorPanel":openColorPanel});
        }
        if (track.type=="bed"){
            return BedTrack({"track":track,"set_value":set_value,"openColorPanel":openColorPanel});
        }
        if (track.type=="bigwig"){
            return BigWigTrack({"track":track,"set_value":set_value,"openColorPanel":openColorPanel});
        }
        if (track.type=="coverage"){
            return CoverageTrack({"track":track,"set_value":set_value,"openColorPanel":openColorPanel});
        }
        if (track.type=="alignments"){
            return AlignmentsTrack({"track":track,"set_value":set_value,"openColorPanel":openColorPanel});
        }
        if (track.type=="basemod_freq"){
            return BasemodfreqTrack({"track":track,"set_value":set_value,"openColorPanel":openColorPanel});
        }
        if (track.type=="hic"){
            return HicTrack({"track":track,"set_value":set_value});
        }
        if (track.type=="sv"){
            return SvTrack({"track":track,"set_value":set_value,"openColorPanel":openColorPanel});
        }
        if (track.type=="copynumber"){
            return CopynumberTrack({"track":track,"set_value":set_value,"openColorPanel":openColorPanel});
        }
    }
    function handleTrackIcon(){
        openTrackTypePanel((t)=>{set_value("type",t)})
        
    }
    return (
        <DndItem id ={track.id} copy_item={()=>copy_track(track.id)} delete_item={()=>delete_track(track.id)} className={className}>
            <div className="trackHeader">
                
            <TrackIcon track_type={track.type} onClick={handleTrackIcon}/>
                
            </div>

            <div className="trackOption">

                

                {track.type!=""?(
                <div className="optionGroup">
                    <div className='formItem'>
                    <label htmlFor={"height"+track.id}>Height (mm):</label>
                    <input id={"height"+track.id} value ={track.height} style={{width:"3em"}} onChange={(e) => set_value("height",e.target.value)}></input>
                    </div>

                    <div className='formItem'>
                        <label htmlFor={"margin_above"+track.id}>Margin above:</label>
                        <input id={"margin_above"+track.id} value ={track.margin_above} style={{width:"2em"}} onChange={(e) => set_value("margin_above",e.target.value)}></input>
                    </div>

                    <div className='formItem'>
                        <label htmlFor={"bounding_box"+track.id}>Box:</label>
                        <input type="checkbox" id={"bounding_box"+track.id} checked={track.bounding_box} onChange={() => set_value("bounding_box",!track.bounding_box)}/>
                    </div>

                    <div className='formItem'>
                        <label htmlFor={"fontscale"+track.id}>Fontscale:</label>
                        <input id={"fontscale"+track.id} value ={track.fontscale} style={{width:"2em"}} onChange={(e) => set_value("fontscale",e.target.value)}></input>
                    </div>

                    <div className='formItem'>
                        <label htmlFor={"label"+track.id}>Label:</label>
                        <input id={"label"+track.id} style={{width:"10em"}} value={track.label}  onChange={(e) => set_value("label",e.target.value)} ></input>
                    </div>
                    <div className='formItem'>
                        <label htmlFor={"label_rotate"+track.id}>Rotate label:</label>
                        <input type="checkbox" id={"label_rotate"+track.id} checked={track.label_rotate} onChange={() => set_value("label_rotate",!track.label_rotate)} ></input>
                    </div>
                </div>

                ):""}

                
                
                {optionTrack()}
                
            </div>

        </DndItem>

    )
}


function ChrTrack({track,set_value}){
    return (
        <>
        <div className="optionGroup">
        <div className='formItem'>
                <label htmlFor={"style"+track.id}>Style:</label>
                <select id={"style"+track.id} value={track.style} onChange={(e) =>{set_value("style",e.target.value)}}> 
                            <option className="dropDownOption" key="default" value="default">default</option>
                            <option className="dropDownOption" key="arrow" value="arrow">arrow</option>
                            <option className="dropDownOption" key="ideogram"  value="ideogram">ideogram</option>
                </select>
            </div>

            <div className='formItem'>
                <label htmlFor={"unit"+track.id}>Unit:</label>
                <select id={"unit"+track.id} value={track.unit} onChange={(e) =>{set_value("unit",e.target.value)}}> 
                            <option className="dropDownOption" key="bp" value="bp">bp</option>
                            <option className="dropDownOption" key="kb" value="kb">kb</option>
                            <option className="dropDownOption" key="Mb"  value="Mb">Mb</option>
                </select>
            </div>

            <div className='formItem'>
                <label htmlFor={"ticklabels_pos"+track.id}>Tick labels position:</label>
                <select id={"ticklabels_pos"+track.id} value={track.ticklabels_pos} onChange={(e) =>{set_value("ticklabels_pos",e.target.value)}}> 
                            <option className="dropDownOption" key="below" value="below">below</option>
                            <option className="dropDownOption" key="above" value="above">above</option>
                            <option className="dropDownOption" key="none"  value="none">none</option>
                </select>
            </div>
            <div className='formItem'>
                <label htmlFor={"ticks_interval"+track.id}>Ticks interval (bp):</label>
                <input id={"ticks_interval"+track.id} style={{width:"4em"}} value={track.ticks_interval}  onChange={(e)=>set_value("ticks_interval",e.target.value)} ></input>
            </div>
        </div>
        
        </>
    )
}

function GenesTrack({track,set_value,openColorPanel}){
    return (
        <>
        <div className="optionGroup">
            <div className='formItem'>
                <label htmlFor={"style"+track.id}>Style:</label>
                <select id={"style"+track.id} value={track.style} onChange={(e) =>{set_value("style",e.target.value)}}> 
                            <option className="dropDownOption" key="default" value="default">default</option>
                            <option className="dropDownOption" key="TSS_arrow" value="TSS_arrow">TSS_arrow</option>
                </select>
            </div>

            <div className='formItem'>
                <label htmlFor={"collapsed"+track.id}>Collapsed:</label>
                <input type="checkbox" id={"collapsed"+track.id} checked={track.collapsed} onChange={() => set_value("collapsed",!track.collapsed)} ></input>
            </div>

            <div className='formItem'>
                <label htmlFor={"only_protein_coding"+track.id}>Only protein coding:</label>
                <input type="checkbox" id={"only_protein_coding"+track.id} checked={track.only_protein_coding} onChange={() => set_value("only_protein_coding",!track.only_protein_coding)} ></input>
            </div>

            <div className='formItem' >
                <label htmlFor={"exon_color"+track.id}>Exon color:</label>
                <ColorButton id={"exon_color"+track.id} color={track.exon_color} setColor={(c)=>set_value("exon_color",c)} openColorPanel={openColorPanel}/>
            </div>

            <div className='formItem'>
                <label htmlFor={"genes"+track.id}>Genes:</label>
                <input id={"genes"+track.id} value={track.genes}  onChange={(e) => set_value("genes",e.target.value)} ></input>
            </div>
        </div>

        </>
    )
}


function BedTrack({track,set_value,openColorPanel}){
    const fileClass= ((!track.file.endsWith(".bed")) && (!track.file.endsWith(".tsv")))? "unsetPath":""; 
    return (
        <>
        <div className="optionGroup">
            <PathEntry id={"file"+track.id} label="File:" value={track.file} set_value={(val) => set_value("file",val)} className={fileClass}/>
            <div className='formItem' >
                <label htmlFor={"color"+track.id} >Color:</label>
                <ColorButton id={"color"+track.id} color={track.color} setColor={(c)=>set_value("color",c)} openColorPanel={openColorPanel}></ColorButton>
            </div>

        </div>
        </>
    )
}


function BigWigTrack({track,set_value,openColorPanel}){
    const fileClass= (track.file.endsWith(".bw") || track.file.endsWith(".bigwig") || track.file.endsWith(".bigWig"))? "":"unsetPath"; 
    return (
        <>
        
            <div className="optionGroup">
                <PathEntry id={"file"+track.id} label="File:" value={track.file} set_value={(val) => set_value("file",val)} className={fileClass} />
                <div className='formItem' >
                    <label htmlFor={"color"+track.id}>Color:</label>
                    <ColorButton id={"color"+track.id} color={track.color} setColor={(c)=>set_value("color",c)} openColorPanel={openColorPanel}></ColorButton>
                </div>
                <div className='formItem'>
                    <label className="withTooltip" htmlFor={"nbins"+track.id}>n_bins:</label>
                    <input id={"nbins"+track.id} value={track.n_bins}  onChange={(e) => set_value("n_bins",e.target.value)} ></input>
                </div>

                {/* scale */}
                <div className='formItem'>
                    <label htmlFor={"scale"+track.id}>Scale:</label>
                    <select id={"scale"+track.id} value={track.scale} onChange={(e) =>{set_value("scale",e.target.value)}}> 
                                <option className="dropDownOption" key="auto" value="auto">auto</option>
                                <option className="dropDownOption" key="auto per region" value="auto per region">auto per region</option>
                                <option className="dropDownOption" key="custom"  value="custom">custom</option>
                    </select>
                </div>
                {(track.scale=="custom")?(
                    <div className='formItem'>
                    <label htmlFor={"scale_max"+track.id}>Max:</label>
                    <input id={"scale_max"+track.id} style={{width:"3em"}} value={track.scale_max} onChange={(e) => set_value("scale_max",e.target.value)} ></input>
                </div>):""
                }
                <div className='formItem'>
                    <label htmlFor={"scale_pos"+track.id}>Scale pos:</label>
                    <select id={"scale_pos"+track.id} value={track.scale_pos} onChange={(e) =>{set_value("scale_pos",e.target.value)}}> 
                                <option className="dropDownOption" key="left" value="left">left</option>
                                <option className="dropDownOption" key="corner" value="corner">corner</option>
                                <option className="dropDownOption" key="corner all"  value="corner all">corner all</option>
                                <option className="dropDownOption" key="none"  value="none">none</option>
                    </select>
                </div>

            </div>
           
        </>
    )
}

function CoverageTrack({track,set_value,openColorPanel}){
    const fileClass= (!track.file.endsWith(".bam"))? "unsetPath":""; 
    return (
        <>
            <div className="optionGroup">
            <PathEntry id={"file"+track.id} label="File:" value={track.file} set_value={(val) => set_value("file",val)} className={fileClass}/>
                <div className='formItem' >
                    <label htmlFor={"color"+track.id}>Color:</label>
                    <ColorButton id={"color"+track.id} color={track.color} setColor={(c)=>set_value("color",c)} openColorPanel={openColorPanel}></ColorButton>
                </div>
                <div className='formItem'>
                    <label className="withTooltip" htmlFor={"nbins"+track.id}>n_bins:</label>
                    <input id={"nbins"+track.id} value={track.n_bins}  onChange={(e) => set_value("n_bins",e.target.value)} ></input>
                </div>

                <div className='formItem'>
                    <label htmlFor={"scale"+track.id}>Scale:</label>
                    <select id={"scale"+track.id} value={track.scale} onChange={(e) =>{set_value("scale",e.target.value)}}> 
                                <option className="dropDownOption" key="auto" value="auto">auto</option>
                                <option className="dropDownOption" key="auto per region" value="auto per region">auto per region</option>
                                <option className="dropDownOption" key="custom"  value="custom">custom</option>
                    </select>
                </div>
                {(track.scale=="custom")?(
                    <div className='formItem'>
                    <label htmlFor={"scale_max"+track.id}>Max:</label>
                    <input id={"scale_max"+track.id} style={{width:"3em"}} value={track.scale_max} onChange={(e) => set_value("scale_max",e.target.value)} ></input>
                </div>):""
                }
                <div className='formItem'>
                    <label htmlFor={"scale_pos"+track.id}>Scale pos:</label>
                    <select id={"scale_pos"+track.id} value={track.scale_pos} onChange={(e) =>{set_value("scale_pos",e.target.value)}}> 
                                <option className="dropDownOption" key="left" value="left">left</option>
                                <option className="dropDownOption" key="corner" value="corner">corner</option>
                                <option className="dropDownOption" key="corner all"  value="corner all">corner all</option>
                                <option className="dropDownOption" key="none"  value="none">none</option>
                    </select>
                </div>

            </div>
            
        </>
    )
}

function AlignmentsTrack({track,set_value,openColorPanel}){
    const fileClass= (!track.file.endsWith(".bam"))? "unsetPath":""; 

    function handleGroupbyChange(e){
        if (e.target_value!=track.group_by){
            set_value("group_by",e.target.value);
            if (e.target.value=="haplotype"){
                set_value("show_unphased",true)
            }
        }
    }

    function setHaplotypeLabel(e,index){
        const Hlabels = [...track.haplotype_labels]
        Hlabels[index] = e.target.value
        set_value("haplotype_labels",Hlabels)
    }
    function setHaplotypeColors(c,index){
        const Hcolors = [...track.haplotype_colors]
        Hcolors[index] = c
        set_value("haplotype_colors",Hcolors)
    }

    function setBasemod(basemod_index,attribute_index,value){
        const basemod= [...track.basemods[basemod_index]];
        basemod[attribute_index] = value;
        const basemods = [...track.basemods];
        basemods[basemod_index] = basemod;
        set_value("basemods",basemods)
    }
    function addBasemod(){
        set_value("basemods",[...track.basemods,["C","h","#ffa500"]])
    }
    function deleteBasemod(){
        set_value("basemods",[track.basemods[0]])
    }
    return (
        <>
        <div className="optionGroup">
            <PathEntry id={"file"+track.id} label="File:" value={track.file} set_value={(val) => set_value("file",val)} className={fileClass} />
            <div className='formItem'>
                <label htmlFor={"hgap_bp"+track.id}>h-gap (bp):</label>
                <input id={"hgap_bp"+track.id} style={{width:"2em"}} value={track.hgap_bp}  onChange={(e)=>set_value("hgap_bp",e.target.value)} ></input>
            </div>

            <div className='formItem'>
                <label htmlFor={"vgap_frac"+track.id}>v-gap (frac):</label>
                <input id={"vgap_frac"+track.id} style={{width:"2em"}} value={track.vgap_frac}  onChange={(e)=>set_value("vgap_frac",e.target.value)} ></input>
            </div>


            <div className='formItem' >
                <label htmlFor={"read_color"+track.id}>Read color:</label>
                <ColorButton id={"read_color"+track.id} color={track.read_color} setColor={(c)=>set_value("read_color",c)} openColorPanel={openColorPanel}/>
            </div>
        </div>

        {/* Splitreads */}
        <div className="optionGroup">

            <div className='formItem'>
                <label htmlFor={"link_splitreads"+track.id}>Link splitreads:</label>
                <input type="checkbox" id={"link_splitreads"+track.id} checked={track.link_splitreads} onChange={() => set_value("link_splitreads",!track.link_splitreads)} ></input>
            </div>

            <div className='formItem' >
                <label htmlFor={"splitread_color"+track.id}>Splitread color:</label>
                <ColorButton id={"splitread_color"+track.id} color={track.splitread_color} setColor={(c)=>set_value("splitread_color",c)} openColorPanel={openColorPanel}/>
            </div>

            <div className='formItem'>
                <label className="withTooltip" htmlFor={"min_splitreads_breakpoints"+track.id}>Min SR: <span className="tooltip">Minimum number of split reads spanning a breakpoint for this breakpoint to be shown.</span></label>
                <input id={"min_splitreads_breakpoints"+track.id} style={{width:"1em"}} value={track.min_splitreads_breakpoints}  onChange={(e)=>set_value("min_splitreads_breakpoints",e.target.value)} />
            </div>
        </div>
            
        {/* Group by*/}
        <div className="optionGroup">
            <div className='formItem'>
                <label htmlFor={"group_by"+track.id}>Group by:</label>
                <select id={"group_by"+track.id} value={track.group_by} onChange={handleGroupbyChange}> 
                            <option className="dropDownOption" key="none" value="none">none</option>
                            <option className="dropDownOption" key="haplotype" value="haplotype">haplotype</option>
                </select>
            </div>
            {(track.group_by=="haplotype")?(
                <>
                <div className='formItem'>
                    <label htmlFor={"show_unphased"+track.id}>Show unphased:</label>
                    <input type="checkbox" id={"show_unphased"+track.id} checked={track.show_unphased} onChange={() => set_value("show_unphased",!track.show_unphased)} ></input>
                </div>

                <div className='formItem'>
                    <label htmlFor={"exchange_haplotypes"+track.id}>Exchange haplotypes:</label>
                    <input type="checkbox" id={"exchange_haplotypes"+track.id} checked={track.exchange_haplotypes} onChange={() => set_value("exchange_haplotypes",!track.exchange_haplotypes)} ></input>
                </div>

                <div className='formItem'>
                    <label htmlFor={"show_haplotype_colors"+track.id}>Show haplotype colors:</label>
                    <input type="checkbox" id={"show_haplotype_colors"+track.id} checked={track.show_haplotype_colors} onChange={() => set_value("show_haplotype_colors",!track.show_haplotype_colors)} ></input>
                </div>

                <div className='formItem'>
                    <label htmlFor={"labelHaplo1"+track.id}>Haplotype 1:</label>
                    <input id={"labelHaplo1"+track.id} style={{width:"6em"}} value={track.haplotype_labels[0]}  onChange={(e)=>setHaplotypeLabel(e,0)} ></input>
                        {track.show_haplotype_colors?(<ColorButton id={"colorHaplo1"+track.id} color={track.haplotype_colors[0]} setColor={(c)=>setHaplotypeColors(c,0)} openColorPanel={openColorPanel}></ColorButton>
                        ):""}
                </div>
                <div className='formItem'>
                    <label htmlFor={"labelHaplo2"+track.id}>Haplotype 2:</label>
                    <input id={"labelHaplo2"+track.id} style={{width:"6em"}} value={track.haplotype_labels[1]}  onChange={(e)=>setHaplotypeLabel(e,1)} ></input>
                        {track.show_haplotype_colors?(<ColorButton id={"colorHaplo2"+track.id} color={track.haplotype_colors[1]} setColor={(c)=>setHaplotypeColors(c,1)} openColorPanel={openColorPanel}></ColorButton>
                        ):""}
                </div>

                <div className='formItem'>
                    {track.show_unphased?(
                    <>
                    <label htmlFor={"labelHaplo3"+track.id}>Unphased:</label>
                    <input id={"labelHaplo3"+track.id} style={{width:"6em"}} value={track.haplotype_labels[2]}  onChange={(e)=>setHaplotypeLabel(e,2)} />
                    </>):""}
                    {(track.show_haplotype_colors&&track.show_unphased)?(
                    <ColorButton id={"colorHaplo3"+track.id} color={track.haplotype_colors[2]} setColor={(c)=>setHaplotypeColors(c,2)} openColorPanel={openColorPanel}></ColorButton>
                        ):""}
                </div>

                </>):""}
            </div>

        {/* Color by*/}
        <div className='optionGroup'>
            <div className='formItem'>
                <label htmlFor={"color_by"+track.id}>Color by:</label>
                <select id={"color_by"+track.id} value={track.color_by} onChange={(e)=>set_value("color_by",e.target.value)}> 
                            <option className="dropDownOption" key="none" value="none">none</option>
                            <option className="dropDownOption" key="basemod" value="basemod">basemod</option>
                </select>
            </div>
            {(track.color_by=="basemod")?(
                <>
                <div className='formItem'>
                    <label htmlFor={"color_unmodified"+track.id}>Unmodified:</label>
                    <ColorButton id={"color_unmodified"+track.id} color={track.color_unmodified} setColor={(c)=>set_value("color_unmodified",c)} openColorPanel={openColorPanel}/>
                </div>

                <div className='formItem'>
                    <label htmlFor={"base1"+track.id}>Basemod 1:</label>
                    <input id={"base1"+track.id} style={{width:"1em"}} value={track.basemods[0][0]}  onChange={(e)=>setBasemod(0,0,e.target.value)}/>
                    <input id={"mod1"+track.id} style={{width:"1em"}} value={track.basemods[0][1]}  onChange={(e)=>setBasemod(0,1,e.target.value)}/>
                    <ColorButton id={"color_basemod1"+track.id} color={track.basemods[0][2]} setColor={(c)=>setBasemod(0,2,c)} openColorPanel={openColorPanel}/>
                </div>
                {track.basemods.length==1? (<button onClick={addBasemod}>Add basemod</button>):""}
                {track.basemods.length==2? (
                <>
                <div className='formItem'>
                    <label htmlFor={"base2"+track.id}>Basemod 2:</label>
                    <input id={"base2"+track.id} style={{width:"1em"}} value={track.basemods[1][0]}  onChange={(e)=>setBasemod(1,0,e.target.value)}/>
                    <input id={"mod2"+track.id} style={{width:"1em"}} value={track.basemods[1][1]}  onChange={(e)=>setBasemod(1,1,e.target.value)}/>
                    <ColorButton id={"color_basemod2"+track.id} color={track.basemods[1][2]} setColor={(c)=>setBasemod(1,2,c)} openColorPanel={openColorPanel}/>
                </div>
                <button onClick={deleteBasemod}>Delete basemod</button>
                
                </>
                ):""}
                <div className='formItem'>
                    <label htmlFor={"fix_hardclip_basemod"+track.id}>Fix hardclip basemod:</label>
                    <input type="checkbox" id={"fix_hardclip_basemod"+track.id} checked={track.fix_hardclip_basemod} onChange={() => set_value("fix_hardclip_basemod",!track.fix_hardclip_basemod)} ></input>
                </div>

                </>
                ):""
            }
            </div>

        </>
    )
}


function HicTrack({track,set_value}){
    const fileClass= ((!track.file.endsWith(".cool")) && (!track.file.includes(".mcool::resolutions/")))? "unsetPath":""; 

    return (
        <>
        <div className="optionGroup">
            <PathEntry id={"file"+track.id} label="File:" value={track.file} set_value={(val) => set_value("file",val)} className={fileClass} />
            <div className='formItem'>
                <label htmlFor={"color_map"+track.id}>Color map:</label>
                <select id={"color_map"+track.id} value={track.color_map} onChange={(e) =>{set_value("color_map",e.target.value)}}> 
                            <option className="dropDownOption" key="red" value="red">red</option>
                            <option className="dropDownOption" key="heat" value="heat">heat</option>
                </select>
            </div>

            <div className='formItem'>
                <label htmlFor={"scale"+track.id}>Scale:</label>
                <select id={"scale"+track.id} value={track.scale} onChange={(e) =>{set_value("scale",e.target.value)}}> 
                            <option className="dropDownOption" key="auto" value="auto">auto</option>
                            <option className="dropDownOption" key="custom"  value="custom">custom</option>
                </select>
            </div>
                {(track.scale=="custom")?(
                    <>
                    <div className='formItem'>
                    <label htmlFor={"scale_min"+track.id}>Min:</label>
                    <input id={"scale_min"+track.id} style={{width:"3em"}} value={track.scale_min} onChange={(e) => set_value("scale_min",e.target.value)} ></input>
                    </div>
                    <div className='formItem'>
                    <label htmlFor={"scale_max"+track.id}>Max:</label>
                    <input id={"scale_max"+track.id} style={{width:"3em"}} value={track.scale_max} onChange={(e) => set_value("scale_max",e.target.value)} ></input>
                    </div>
                    </>):(
                        <div className='formItem'>
                        <label htmlFor={"scale_max_percentile"+track.id}>Max percentile:</label>
                        <input id={"scale_max_percentile"+track.id} style={{width:"3em"}} value={track.scale_max_percentile} onChange={(e) => set_value("scale_max_percentile",e.target.value)} ></input>
                        </div>
                    )
                }
            <div className='formItem'>
                <label htmlFor={"show_colorbar"+track.id}>Show colorbar:</label>
                <input type="checkbox" id={"show_colorbar"+track.id} checked={track.show_colorbar} onChange={() => set_value("show_colorbar",!track.show_colorbar)} />
            </div>
            
        </div>

        <div className="optionGroup">

            <div className='formItem'>
                <label htmlFor={"max_dist"+track.id}>Max dist (kb):</label>
                <input id={"max_dist"+track.id} style={{width:"4em"}} value={track.max_dist}  onChange={(e) => set_value("max_dist",e.target.value)} />
            </div>

            <div className='formItem'>
                <label htmlFor={"extend"+track.id}>Extend:</label>
                <input type="checkbox" id={"extend"+track.id} checked={track.extend} onChange={() => set_value("extend",!track.extend)} />
            </div>
            <div className='formItem'>
                <label htmlFor={"interactions_across_regions"+track.id}>Interactions across regions:</label>
                <input type="checkbox" id={"interactions_across_regions"+track.id} checked={track.interactions_across_regions} onChange={() => set_value("interactions_across_regions",!track.interactions_across_regions)} />
            </div>
            {track.interactions_across_regions?(
                <div className='formItem'>
                <label htmlFor={"double_interactions_across_regions"+track.id}>Double across regions:</label>
                <input type="checkbox" id={"double_interactions_across_regions"+track.id} checked={track.double_interactions_across_regions} onChange={() => set_value("double_interactions_across_regions",!track.double_interactions_across_regions)} />
            </div>
            ):""}
            <div className='formItem'>
                <label htmlFor={"pixel_border"+track.id}>Pixel border:</label>
                <input type="checkbox" id={"pixel_border"+track.id} checked={track.pixel_border} onChange={() => set_value("pixel_border",!track.pixel_border)} />
            </div>
            <div className='formItem'>
                <label htmlFor={"upside_down"+track.id}>Upside down:</label>
                <input type="checkbox" id={"upside_down"+track.id} checked={track.upside_down} onChange={() => set_value("upside_down",!track.upside_down)} />
            </div>
        </div>
            
        </>
    )
}

function SvTrack({track,set_value,openColorPanel}){
    const fileClass= (track.file.endsWith(".vcf") || track.file.endsWith(".vcf.gz") ||track.file.endsWith(".tsv"))? "":"unsetPath"; 
    return (
        <>
        <div className="optionGroup">
            <PathEntry id={"file"+track.id} label="File:" value={track.file} set_value={(val) => set_value("file",val)} className={fileClass} />

            <div className='formItem'>
                <label htmlFor={"lw"+track.id}>Line width:</label>
                <input id={"lw"+track.id} style={{width:"2em"}} value={track.lw}  onChange={(e) => set_value("lw",e.target.value)} ></input>
            </div>

            <div className='formItem'>
                <label htmlFor={"color_del"+track.id}>DEL-like:</label>
                <ColorButton id={"color_del"+track.id} color={track.color_del} setColor={(c)=>set_value("color_del",c)} openColorPanel={openColorPanel}/>
            </div>
            <div className='formItem'>
                <label htmlFor={"color_dup"+track.id}>DUP-like:</label>
                <ColorButton id={"color_dup"+track.id} color={track.color_dup} setColor={(c)=>set_value("color_dup",c)} openColorPanel={openColorPanel}/>
            </div>
            <div className='formItem'>
                <label htmlFor={"color_h2h"+track.id}>H2H:</label>
                <ColorButton id={"color_h2h"+track.id} color={track.color_h2h} setColor={(c)=>set_value("color_h2h",c)} openColorPanel={openColorPanel}/>
            </div>
            <div className='formItem'>
                <label htmlFor={"color_t2t"+track.id}>T2T:</label>
                <ColorButton id={"color_t2t"+track.id} color={track.color_t2t} setColor={(c)=>set_value("color_t2t",c)} openColorPanel={openColorPanel}/>
            </div>
            <div className='formItem'>
                <label htmlFor={"color_trans"+track.id}>Translocation:</label>
                <ColorButton id={"color_trans"+track.id} color={track.color_trans} setColor={(c)=>set_value("color_trans",c)} openColorPanel={openColorPanel}/>
            </div>


        </div>
            
            
            

        </>
    )
}


function CopynumberTrack({track,set_value,openColorPanel}){
    const fileClass= (track.freec_ratios!=="" || track.freec_CNAs!=="" || track.purple_cn!=="")? "":"unsetPath"; 
    return (
        <>
        <div className="optionGroup">
            <PathEntry id={"freec_ratios"+track.id} label="Ratios file:" value={track.freec_ratios} set_value={(val) => set_value("freec_ratios",val)} className={fileClass}/>
            <PathEntry id={"freec_CNAs"+track.id} label="CNAs file:" value={track.freec_CNAs} set_value={(val) => set_value("freec_CNAs",val)} className={fileClass}/>
            <PathEntry id={"purple_cn"+track.id} label="Purple CN file:" value={track.purple_cn} set_value={(val) => set_value("purple_cn",val)} className={fileClass}/>

            <div className='formItem'>
                <label htmlFor={"genes"+track.id}>Genes:</label>
                <input id={"genes"+track.id} value={track.genes}  onChange={(e) => set_value("genes",e.target.value)} />
            </div>

            <div className='formItem'>
                <label htmlFor={"min_cn"+track.id}>Min CN:</label>
                <input id={"min_cn"+track.id} style={{width:"2em"}} value={track.min_cn}  onChange={(e) => set_value("min_cn",e.target.value)}/>
            </div>

            <div className='formItem'>
                <label htmlFor={"max_cn"+track.id}>Max CN:</label>
                <input id={"max_cn"+track.id} style={{width:"2em"}} value={track.max_cn}  onChange={(e) => set_value("max_cn",e.target.value)}/>
            </div>
        </div>

        <div className="optionGroup">
            {/* Colors */}
            <div className='formItem'>
                <label htmlFor={"color_normal"+track.id}>Normal:</label>
                <ColorButton id={"color_normal"+track.id} color={track.color_normal} setColor={(c)=>set_value("color_normal",c)} openColorPanel={openColorPanel}/>
            </div>
            <div className='formItem'>
                <label htmlFor={"color_loss"+track.id}>Loss:</label>
                <ColorButton id={"color_loss"+track.id} color={track.color_loss} setColor={(c)=>set_value("color_loss",c)} openColorPanel={openColorPanel}/>
            </div>
            <div className='formItem'>
                <label htmlFor={"color_gain"+track.id}>Gain:</label>
                <ColorButton id={"color_gain"+track.id} color={track.color_gain} setColor={(c)=>set_value("color_gain",c)} openColorPanel={openColorPanel}/>
            </div>
            <div className='formItem'>
                <label htmlFor={"color_cnloh"+track.id}>CNLOH:</label>
                <ColorButton id={"color_cnloh"+track.id} color={track.color_cnloh} setColor={(c)=>set_value("color_cnloh",c)} openColorPanel={openColorPanel}/>
            </div>
        </div>

        <div className="optionGroup">

            {/* Grid */}
            <div className='formItem'>
                <label htmlFor={"grid"+track.id}>Grid:</label>
                <input type="checkbox" id={"grid"+track.id} checked={track.grid} onChange={() => set_value("grid",!track.grid)} ></input>
            </div>
            {track.grid? (
                <>
                <div className='formItem'>
                    <label htmlFor={"grid_major"+track.id}>Grid major:</label>
                    <input type="checkbox" id={"grid_major"+track.id} checked={track.grid_major} onChange={() => set_value("grid_major",!track.grid_major)} ></input>
                </div>
                <div className='formItem'>
                    <label htmlFor={"grid_minor"+track.id}>Grid minor:</label>
                    <input type="checkbox" id={"grid_minor"+track.id} checked={track.grid_minor} onChange={() => set_value("grid_minor",!track.grid_minor)} ></input>
                </div>
                <div className='formItem'>
                    <label htmlFor={"grid_cn"+track.id}>Grid CN:</label>
                    <input type="checkbox" id={"grid_cn"+track.id} checked={track.grid_cn} onChange={() => set_value("grid_cn",!track.grid_cn)} ></input>
                </div>
            </>
            ):""}

        </div>
            

        </>
    )
}


export function TrackIcon({track_type, onClick}){
    return(
        <div className="trackIcon" onClick={onClick}>

        <div style={{backgroundColor:"white",borderRadius:"30px",display:"flex",alignItems:"center",justifyContent:"center",alignContent:"center"}}>
        {(track_type=="chr_axis")?(
            <svg width="15mm" height="15mm" version="1.1" viewBox="0 0 10 10" xmlns="http://www.w3.org/2000/svg">
            <g>
             <rect x="1" y="4.8" width="8" height=".4" opacity=".97208"/>
             <rect x="1" y="4.35" width=".2" height="1.3" opacity=".97208"/>
             <rect x="8.98" y="4.35" width=".2" height="1.3" opacity=".97208"/>
             <rect x="3" y="4.35" width=".2" height="1.3" opacity=".97208"/>
             <rect x="5" y="4.35" width=".2" height="1.3" opacity=".97208"/>
             <rect x="7" y="4.35" width=".2" height="1.3" opacity=".97208"/>
            </g>
           </svg>
        ):""}

        {(track_type=="genes")?(
            <svg width="15mm" height="15mm" version="1.1" viewBox="0 0 10 10" xmlns="http://www.w3.org/2000/svg">
            <g strokeLinecap="round">
             <rect x="1" y="4.75" width="7" height=".5" fill="#ccc" opacity=".97208" strokeWidth=".26834"/>
             <g fill="#2980b9">
              <rect x="1" y="3.85" width=".8" height="2.3" opacity=".97208" strokeWidth=".60204"/>
              <rect x="3" y="3.85" width=".3" height="2.3" opacity=".97208" strokeWidth=".36868"/>
              <rect transform="scale(-1,1)" x="-5.95" y="3.85" width="2" height="2.3" opacity=".97208" strokeWidth=".95192"/>
              <path d="m7.7415 3.85h0.34323l0.91527 1.15-0.91527 1.15h-0.34323z" opacity=".97208" strokeWidth=".507"/>
             </g>
            </g>
           </svg>
        ):""}

        {(track_type=="bed")?(
            <svg width="15mm" height="15mm" version="1.1" viewBox="0 0 10 10" xmlns="http://www.w3.org/2000/svg">
            <g fill="#2980b9" strokeLinecap="round">
             <rect transform="scale(-1,1)" x="-3.1605" y="4.25" width="1.6605" height="1.5" opacity=".97208" strokeWidth="1.3665"/>
             <rect transform="scale(-1,1)" x="-6.1252" y="4.25" width="1.887" height="1.5" opacity=".97208" strokeWidth="1.4567"/>
             <rect transform="scale(-1,1)" x="-8.5" y="4.25" width="1.1101" height="1.5" opacity=".97208" strokeWidth="1.1173"/>
            </g>
           </svg>
        ):""}
        
        {(track_type=="bigwig")?(
            <svg width="15mm" height="15mm" version="1.1" viewBox="0 0 10 10" xmlns="http://www.w3.org/2000/svg">
                <path d="m9 8.4708h-8v-0.1l0.81733-0.26117 0.5365-4.0946 0.92877-2.5327 0.62713 2.2596 0.52704-0.4133 0.40616-0.94285 0.38526 1.6655 0.64971 3.6571 0.28678-0.018519 0.33533-0.96747 0.96653 1.5622 1.5335 0.086247z" fill="#2980b9" opacity=".97208"/>
            </svg>
        ):""}

        {(track_type=="coverage")?(
            <svg width="15mm" height="15mm" version="1.1" viewBox="0 -1 10 11" xmlns="http://www.w3.org/2000/svg">
            <path d="m8.5 8.5h-7v-7h0.7v0.7h0.7v-0.35h0.7v4.025h0.7v-1.05h0.7v0.7h0.7v-1.1375h0.7v-2.5375h0.7v-0.2625h0.70002v0.33875l0.69998-1.275e-4z" fill="#888" opacity=".97208" strokeLinecap="round" strokeWidth=".42835"/>
           </svg>
        ):""}
        {(track_type=="alignments")?(
            <svg width="15mm" height="15mm" version="1.1" viewBox="0 0 10 10" xmlns="http://www.w3.org/2000/svg">
            <g fill="#888" strokeLinecap="round" strokeWidth=".6145">
             <path d="m1.7028 6.7407h4.9264l-2.086e-4 -0.38892 0.88834 0.64821-0.88813 0.64821v-0.38892h-4.9264z" opacity=".97208"/>
             <path d="m1.8293 2.7407h4.9264l-2.086e-4 -0.38892 0.88834 0.64821-0.88813 0.64821v-0.38892h-4.9264z" opacity=".97208"/>
             <path d="m8.2972 4.7407h-4.9264l2.086e-4 -0.38892-0.88834 0.64821 0.88813 0.64821v-0.38892h4.9264z" opacity=".97208"/>
            </g>
           </svg>
        ):""}
        {(track_type=="hic")?(
            <svg width="15mm" height="15mm" version="1.1" viewBox="0 0 10 10" xmlns="http://www.w3.org/2000/svg">
            <path d="m8.5 7.2488h-7l1.25-4.0789 1.25 3.9289 1-2.9176 1 2.9676 1.25-4.3976z" fill="#e74c3c" opacity=".97208" strokeLinecap="round" strokeWidth=".474"/>
           </svg>
        ):""}
        {(track_type=="copynumber")?(
            <svg width="15mm" height="15mm" version="1.1" viewBox="0 0 11 10" xmlns="http://www.w3.org/2000/svg">
             <g strokeLinecap="round" strokeWidth=".474">
              <circle cx="1.7" cy="5" r=".5" opacity=".97208"/>
              <circle cx="2.7" cy="3.125" r=".5" fill="#e55039" opacity=".97208"/>
              <g>
               <circle cx="3.7" cy="5" r=".5" opacity=".97208"/>
               <circle cx="4.7" cy="5" r=".5" opacity=".97208"/>
               <circle cx="5.7" cy="5" r=".5" opacity=".97208"/>
              </g>
              <g fill="#4a69bd">
               <circle cx="8.7" cy="7.625" r=".5" opacity=".97208"/>
               <circle cx="7.7" cy="7.625" r=".5" opacity=".97208"/>
               <circle cx="6.7" cy="7.625" r=".5" opacity=".97208"/>
              </g>
             </g>
            </svg>
        ):""}

        {(track_type=="sv")?(
            <svg width="15mm" height="15mm" version="1.1" viewBox="0 0 10 10" xmlns="http://www.w3.org/2000/svg">
            <g strokeLinecap="round">
             <path d="m2 6c0.030346-0.04389 0.81423-4.0059 3-4 2.2109 0.00601 3 4 3 4" fill="none" opacity=".97208" stroke="#27ae60" strokeWidth=".5"/>
            </g>
           </svg>
        ):""}

        {(track_type=="basemod_freq")?(
            <svg width="15mm" height="15mm" version="1.1" viewBox="0 0 10 10" xmlns="http://www.w3.org/2000/svg">
            <g fill="none" strokeLinecap="round">
             <path d="m1.034 3.3054c1.2872-0.15163 0.51235 2.4665 1.2218 3.5084 0.23162 0.34177 1.1852 0.15204 1.5823 0.1485 1.1557-0.010299 0.47151-2.9717 1.3575-3.7475 0.41083-0.37614 2.7615-0.18619 3.8657-0.17441" opacity=".97208" stroke="#27ae60" strokeWidth=".5"/>
             <path d="m0.93871 3.1564c1.9901-0.039321 5.146-0.22674 5.547 0.17847 0.42505 0.42951 0.57483 1.2406 1.0938 1.288 0.66505 0.060831 1.1817-1.1639 1.4669-1.5903" opacity=".97208" stroke="#e67e22" strokeWidth=".5"/>
            </g>
           </svg>
        ):""}
        </div>
        


        <div style={{color:"#ffffff", textAlign:"center"}}>
            {track_type}
            {track_type==""? "click to select track type":""}
            </div>
        </div>     
        
    )
    
}


function select_file(set_value,initial_value){
    fetch('/browse',{headers: {'Content-Type': 'application/json'}, body: JSON.stringify({path:initial_value}),
    method:"POST"}).then(res => res.json()).then(data => {
        if (data.path.length>0) set_value(data.path);
      });
}
export function PathEntry({id,label,value,set_value,className=""}){
    return(
        <div className='formItem '>
            <label htmlFor={id}>{label}</label>
            <input id={id} value={value} className={className} onChange={(e)=>set_value(e.target.value) } style={{marginLeft:"2mm", marginRight:"1mm"}}></input>
            <button name="file" onClick={()=>select_file(set_value,value)} className="buttonImage">
            <svg width="4.6mm" height="4.5mm" version="1.1" viewBox="0 0 576 512" xmlns="http://www.w3.org/2000/svg">
                <path d="M384 480h48c11.4 0 21.9-6 27.6-15.9l112-192c5.8-9.9 5.8-22.1 .1-32.1S555.5 224 544 224H144c-11.4 0-21.9 6-27.6 15.9L48 357.1V96c0-8.8 7.2-16 16-16H181.5c4.2 0 8.3 1.7 11.3 4.7l26.5 26.5c21 21 49.5 32.8 79.2 32.8H416c8.8 0 16 7.2 16 16v32h48V160c0-35.3-28.7-64-64-64H298.5c-17 0-33.3-6.7-45.3-18.7L226.7 50.7c-12-12-28.3-18.7-45.3-18.7H64C28.7 32 0 60.7 0 96V416c0 35.3 28.7 64 64 64H87.7 384z" fill="#555555"/>
            </svg>
            </button>
        </div>
    )
}