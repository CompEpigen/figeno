import {React, Fragment,useState,useRef, useCallback} from 'react';
import useClickOutside from "./useClickOutside";
import "./style.css";
import { TrackIcon } from './Track';


export function TrackTypePanel({setTrackType,close}){
    const track_types = ["chr_axis", "genes", "bed", "bigwig", "coverage", "alignments", "basemod_freq", "hic", "sv", "copynumber"];
    const popover = useRef();
    useClickOutside(popover, close);

    return(
    <div className="TrackTypePanel" ref={popover}>
        <h2> Select track type </h2>
        <div className="TrackTypeGrid">
            {track_types.map((t)=>{
                return <TrackIcon key={t} track_type={t} onClick={()=>{setTrackType(t);close();}}/>
            })}
        </div>
        <button onClick={close}>Cancel</button>
    </ div>
    )
}
