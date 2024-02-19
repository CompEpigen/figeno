import React from 'react';
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
        unit:"kb",
        ticklabels_pos:"below",
        ticks_interval:"auto",
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
        genes:"auto"
    },
    "bed":{
        height:"10",
        margin_above:"1.5",
        bounding_box:false,
        fontscale:1.0,
        file:"",
        color:"#444444",
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
        scale_pos:"corner"
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
        scale_pos:"corner"
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
        splitread_color:"#999999",
        link_splitreads:false,
        min_splitreads_breakpoints:2,

        group_by:"none",
        show_unphased:true,
        exchange_haplotypes:false,
        show_haplotype_colors:true,
        haplotype_colors:["#27ae60","#e67e22","#808080"],
        haplotype_labels:["Haplotype 1", "Haplotype 2", "Unphased"],

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
        bams:[]
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
        color_trans:"#27ae60"
    },
    "copynumber":{
        height:"30",
        margin_above:"1.5",
        bounding_box:true,
        fontscale:1.0,
        label:"",
        label_rotate:false,
        freec_ratios:"",
        freec_CNAs:"",
        purple_cn:"",
        genes:"",
        min_cn:"",
        max_cn:"",
        grid: true,
        grid_major:true,
        grid_minor:true,
        grid_cn:true,
        color_normal:"#000000",
        color_loss:"#4a69bd",
        color_gain:"#e55039",
        color_cnloh:"#f6b93b"
    },
    "other":{
        height:"8",
        margin_above:"1.5",
        bounding_box:false,
        fontscale:1.0
    }
}

export function TracksContainer({tracksList,setTracksList,openColorPanel, openTrackTypePanel}) {

    

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
        const newTrack = {...tracksList[index],id:crypto.randomUUID()};
        const newList = [...tracksList]
        newList.splice(index,0,newTrack);
        setTracksList(newList);
    }

    function create_track({id,type,track={}}){
        let new_track = {id:id, type:type};
        if (["chr_axis","genes","bed","bigwig","coverage","alignments", "basemod_freq","hic","sv","copynumber"].includes(new_track.type)){
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
        setTracksList([create_track({id:crypto.randomUUID(),type:new_type,track:track}),...tracksList]);
    }
    function getTrackById (tracks, id) {
        return tracks.find((track) => track.id === id);
    }
    function delete_track(id){
        setTracksList((tracks_list) => {return tracks_list.filter(track => track.id!==id)})
    }

    function open_files(){
        fetch('/open_files').then(res => res.json()).then(data => {   
        if (!data.hasOwnProperty("files")){return }
          const tracks=[];
          for (const f of data.files){
            if (f.endsWith(".bed")){
                tracks.push({id:crypto.randomUUID(),type:"bed",...defaultTrackValues["bed"],file:f})
            }
            else if (f.endsWith(".bw") || (f.endsWith(".bigwig") || (f.endsWith("bigWig")))){
                tracks.push({id:crypto.randomUUID(),type:"bigwig",...defaultTrackValues["bigwig"],file:f})
            }
            else if (f.endsWith(".cool")){
                tracks.push({id:crypto.randomUUID(),type:"hic",...defaultTrackValues["hic"],file:f})
            }
            else if (f.endsWith(".bam")){
                tracks.push({id:crypto.randomUUID(),type:"alignments",...defaultTrackValues["alignments"],file:f})
            }
          }
          setTracksList([...tracks,...tracksList]);
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
                    return <Track key={track.id} track={track} set_value={(attribute,value)=>set_value(track.id,attribute,value)}  className={"track"} copy_track={copy_track} delete_track={delete_track} openColorPanel={openColorPanel} openTrackTypePanel={openTrackTypePanel} />
                })}
        </DndContainer>
        </>
    )
}