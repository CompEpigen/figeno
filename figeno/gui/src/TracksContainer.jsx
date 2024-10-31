import React from 'react';
import {v4 as uuid4} from 'uuid';
import "./style.css";
import { DndContainer } from './DndContainer';
import {Track} from './Track'


export const defaultTrackValues={
    "chr_axis":{
        height:"10",
        margin_above:"1.5",
        bounding_box:false,
        fontscale:1.0,
        label:"",
        label_rotate:false,
        style:"default",
        lw_scale:"1.0",
        ticklabels_pos:"below",
        unit:"kb",
        ticks_interval:"auto",
        ticks_angle:0,
        chr_prefix: "chr"
    },
    "genes":{
        height:"10",
        margin_above:"1.5",
        bounding_box:false,
        fontscale:1.0,
        label:"",
        label_rotate:false,
        style:"default",
        collapsed:true,
        only_protein_coding:true,
        exon_color:"#2980b9",
        genes:"auto",
        show_gene_names:true
    },
    "bed":{
        height:"3",
        margin_above:"1.5",
        bounding_box:false,
        fontscale:1.0,
        file:"",
        color:"#444444",
        show_names:false,
        show_strand:false,
        label:"",
        label_rotate:false
    },
    "bigwig":{
        height:"10",
        margin_above:"1.5",
        bounding_box:false,
        fontscale:1.0,
        label:"",
        label_rotate:false,
        file:"",
        color:"#2980b9",
        n_bins:500,
        scale:"auto",
        scale_max:"",
        group:"1",
        scale_pos:"corner",
        upside_down:false
    },
    "coverage":{
        height:"10",
        margin_above:"1.5",
        bounding_box:false,
        fontscale:1.0,
        label:"",
        label_rotate:false,
        file:"",
        color:"#888888",
        n_bins:500,
        scale:"auto",
        scale_max:"",
        group:"1",
        scale_pos:"corner",
        upside_down:false
    },
    "alignments":{
        height:"50",
        margin_above:"1.5",
        bounding_box:false,
        fontscale:1.0,
        file:"",
        label:"",
        label_rotate:false,
        hgap_bp:30,
        vgap_frac:0.3,
        read_color:"#cccccc",
        link_splitreads:false,
        splitread_color:"#999999",
        link_color:"#999999",
        link_lw:0.2,
        only_show_splitreads:false,
        min_splitreads_breakpoints:2,
        

        group_by:"none",
        show_unphased:true,
        exchange_haplotypes:false,
        show_haplotype_colors:true,
        haplotype_colors:["#27ae60","#e67e22","#808080"],
        haplotype_labels:["HP1", "HP2", "Unphased"],

        color_by:"none",
        color_unmodified:"#0f57e5",
        basemods:[["C","m","#f40202"]],
        fix_hardclip_basemod:false
    },
    "basemod_freq":{
        height:"20",
        margin_above:"1.5",
        bounding_box:true,
        fontscale:1.0,
        label:"Methylation freq",
        label_rotate:true,
        style:"lines",
        smooth:4,
        gap_frac:0.1,
        bams:[],
        bedmethyls:[]
    },
    "hic":{
        height:"50",
        margin_above:"1.5",
        bounding_box:true,
        fontscale:1.0,
        label:"",
        label_rotate:false,
        file:"",
        color_map:"red",
        pixel_border:false,
        upside_down:false,
        max_dist:700,
        extend:true,
        scale:"auto",
        scale_max_percentile:90,
        scale_min:"0.0",
        scale_max:"0.15",
        show_colorbar: false,
        interactions_across_regions:true,
        double_interactions_across_regions:true
    },
    "sv":{
        height:"15",
        margin_above:"1.5",
        bounding_box:true,
        fontscale:1.0,
        label:"",
        label_rotate:false,
        file:"",
        lw:"0.5",
        color_del:"#4a69bd",
        color_dup:"#e55039",
        color_t2t:"#8e44ad",
        color_h2h:"#8e44ad",
        color_trans:"#27ae60",
        min_sv_height:0.1
    },
    "copynumber":{
        height:"30",
        margin_above:"1.5",
        bounding_box:true,
        fontscale:1.0,
        label:"",
        label_rotate:false,
        input_type:"freec",
        freec_ratios:"",
        freec_CNAs:"",
        purple_cn:"",
        delly_cn:"",
        delly_CNAs:"",
        genes:"",
        ploidy:"2",
        min_cn:"",
        max_cn:"",
        grid: true,
        grid_major:true,
        grid_minor:true,
        grid_cn:true,
        marker_size:"0.7",
        color_normal:"#000000",
        color_loss:"#4a69bd",
        color_gain:"#e55039",
        color_cnloh:"#f6b93b"
    },
    "ase":{
        height:"50",
        margin_above:"1.5",
        bounding_box:false,
        fontscale:1.0,
        label:"RNA,DNA",
        label_rotate:true,
        file:"",
        vcf_DNA:"",
        min_depth:"6",
        color1:"#e55039",
        color2:"#4a69bd",
        max_bar_width:"10.0",
        lw: "0.1",
        only_exonic: false,
        grid:false
    },
    "other":{
        height:"8",
        margin_above:"1.5",
        bounding_box:false,
        fontscale:1.0
    }
}

export function TracksContainer({tracksList,setTracksList,openColorPanel, openTrackTypePanel, setFileDialogData,setFileDialogActive}) {

    

    function set_value(id,attribute,value){
        setTracksList((tracks)=>{
            return tracks.map((track)=>{
                if (track.id===id) {
                    track[attribute]=value;
                    if (attribute==="type"){
                        track = create_track({id:track.id,type:value})
                    }
                    return track;}
                else {return track;}
            })
        });
        //track[attribute] = value
    }
    
    function get_tracks(){
        return tracksList;
    }
    function copy_track(id){
        const index = tracksList.findIndex((t)=>t.id==id);
        const newTrack = {...tracksList[index],id:uuid4()};
        const newList = [...tracksList]
        newList.splice(index,0,newTrack);
        setTracksList(newList);
    }

    function create_track({id,type,track={}}){
        let new_track = {id:id, type:type};
        if (["chr_axis","genes","bed","bigwig","coverage","alignments", "basemod_freq","hic","sv","copynumber","ase"].includes(new_track.type)){
            for (const attribute in defaultTrackValues[new_track.type]){
                new_track[attribute] = track.hasOwnProperty(attribute) ? track[attribute] : defaultTrackValues[new_track.type][attribute]
            }
        }
        else{
            for (const attribute in defaultTrackValues["other"]){
                new_track[attribute] = track.hasOwnProperty(attribute) ? track[attribute] : defaultTrackValues["other"][attribute]
            }
        }
        return new_track
    }
    
    function handle_add_track(){
        openTrackTypePanel((t)=>{add_track({"type":t})});
    }

    function add_track(track={}){
        const new_type = track.hasOwnProperty("type") ? track.type : "";
        setTracksList([create_track({id:uuid4(),type:new_type,track:track}),...tracksList]);
    }
    function getTrackById (tracks, id) {
        return tracks.find((track) => track.id === id);
    }
    function delete_track(id){
        setTracksList((tracks_list) => {return tracks_list.filter(track => track.id!==id)})
    }

    function add_tracks_files(files){
        const tracks=[];
        for (const f of files){
            if (f.endsWith(".bed")){
                tracks.push({id:uuid4(),type:"bed",...defaultTrackValues["bed"],file:f})
            }
            else if (f.endsWith(".bw") || (f.endsWith(".bigwig") || (f.endsWith("bigWig")))){
                tracks.push({id:uuid4(),type:"bigwig",...defaultTrackValues["bigwig"],file:f})
            }
            else if (f.endsWith(".cool")){
                tracks.push({id:uuid4(),type:"hic",...defaultTrackValues["hic"],file:f})
            }
            else if (f.endsWith(".bam")){
                tracks.push({id:uuid4(),type:"alignments",...defaultTrackValues["alignments"],file:f})
            }
        }
        setTracksList([...tracks,...tracksList]);
    }

    function open_files(){
        fetch('/browse',{headers: {'Content-Type': 'application/json'}, body: JSON.stringify({path:"",dialog_type:"open_files"}),
    method:"POST"}).then(res => res.json()).then(data => {  
            let files=[];
            if (data.hasOwnProperty("current_dir")){
                setFileDialogData({current_dir:data.current_dir,dirs:data.dirs,files:data.files,dialog_type:"open_files",update:function(f){add_tracks_files([f])}})
                setFileDialogActive(true);
            }
            else if (data.hasOwnProperty("files")) add_tracks_files(data.files);
        });
      }


    const header=<div className="TrackContainerHeader">
        <h3>Tracks</h3>
    <button className="ButtonHeader" onClick={(e)=>handle_add_track()}>
    <svg xmlns="http://www.w3.org/2000/svg" height="24" viewBox="0 -960 960 960" width="24"><path d="M440-440H200v-80h240v-240h80v240h240v80H520v240h-80v-240Z" fill="#ffffff" /></svg>
        Add track
    </button>
    <button className="ButtonHeader" onClick={(e)=>open_files()}>
    <svg style={{paddingRight:"10px"}} width="20" height="18" version="1.1" viewBox="0 0 576 512" xmlns="http://www.w3.org/2000/svg">
                <path d="M384 480h48c11.4 0 21.9-6 27.6-15.9l112-192c5.8-9.9 5.8-22.1 .1-32.1S555.5 224 544 224H144c-11.4 0-21.9 6-27.6 15.9L48 357.1V96c0-8.8 7.2-16 16-16H181.5c4.2 0 8.3 1.7 11.3 4.7l26.5 26.5c21 21 49.5 32.8 79.2 32.8H416c8.8 0 16 7.2 16 16v32h48V160c0-35.3-28.7-64-64-64H298.5c-17 0-33.3-6.7-45.3-18.7L226.7 50.7c-12-12-28.3-18.7-45.3-18.7H64C28.7 32 0 60.7 0 96V416c0 35.3 28.7 64 64 64H87.7 384z" fill="#ffffff"/>
            </svg>
        Open files
    </button>
    </div>;

    function show_active_item(active_id){
        return <Track track={{...getTrackById(tracksList,active_id),id:active_id+"overlay"}} className={"track trackOverlay"} /> 
    }

    return (
        <>
        <DndContainer header={header} items={tracksList} setItems={setTracksList} show_active_item={show_active_item}>
        {tracksList.map((track)=>{
                    return <Track key={track.id} track={track} set_value={(attribute,value)=>set_value(track.id,attribute,value)}  className={"track"} copy_track={copy_track} delete_track={delete_track} openColorPanel={openColorPanel} openTrackTypePanel={openTrackTypePanel} setFileDialogActive={setFileDialogActive} setFileDialogData={setFileDialogData} />
                })}
        </DndContainer>
        </>
    )
}