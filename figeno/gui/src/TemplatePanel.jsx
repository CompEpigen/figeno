import {React, Fragment,useState,useRef, useCallback} from 'react';
import useClickOutside from "./useClickOutside";
import "./style.css";
import { TrackIcon } from './Track';
import { defaultTrackValues } from './TracksContainer';


export function TemplatePanel({setTracksList,setRegionsList,setHighlightsList,setGeneralParams,close}){
    const templates = ["bigwig", "hic", "asm", "wgs_chr", "wgs_circos"];
    const popover = useRef();
    useClickOutside(popover, close);

    function select_bigwig(){
        setGeneralParams((params)=> {return {...params,layout:"horizontal"}});
        setRegionsList([{id:crypto.randomUUID(),chr:"",start:"",end:"",color:"#f4a460"}]);
        let bigwigTrack= {id:crypto.randomUUID(),type:"bigwig",...defaultTrackValues.bigwig};
        let genesTrack={id:crypto.randomUUID(),type:"genes",...defaultTrackValues.genes};
        let chrTrack={id:crypto.randomUUID(),type:"chr_axis",...defaultTrackValues["chr_axis"]};
        setHighlightsList([]);
        setTracksList([bigwigTrack,genesTrack,chrTrack])
        close();
    }

    function select_hic(){
        setGeneralParams((params)=> {return {...params,layout:"horizontal"}});
        setRegionsList([{id:crypto.randomUUID(),chr:"",start:"",end:"",color:"#f4a460"}]);
        let hicTrack= {id:crypto.randomUUID(),type:"hic",...defaultTrackValues.hic};
        let genesTrack={id:crypto.randomUUID(),type:"genes",...defaultTrackValues.genes};
        let chrTrack={id:crypto.randomUUID(),type:"chr_axis",...defaultTrackValues["chr_axis"]};
        setHighlightsList([]);
        setTracksList([hicTrack,genesTrack,chrTrack])
        close();
    }

    function select_asm(){
        setGeneralParams((params)=> {return {...params,layout:"horizontal"}});
        setRegionsList([{id:crypto.randomUUID(),chr:"",start:"",end:"",color:"#f4a460"}]);
        let alignmentsTrack= {id:crypto.randomUUID(),type:"alignments",...defaultTrackValues.alignments,group_by:"haplotype",color_by:"basemod"};
        let basemodfreqTrack= {id:crypto.randomUUID(),type:"basemod_freq",...defaultTrackValues["basemod_freq"],
                        bams:[{id:crypto.randomUUID(),"file": "","base": "C","mod": "m","min_coverage": 6,"linewidth": 3,"opacity": 1,
                        "fix_hardclip": false,"split_by_haplotype": true,"colors": ["#27ae60","#e67e22"]}]};
        let genesTrack={id:crypto.randomUUID(),type:"genes",...defaultTrackValues.genes};
        let chrTrack={id:crypto.randomUUID(),type:"chr_axis",...defaultTrackValues["chr_axis"]};
        setHighlightsList([]);
        setTracksList([alignmentsTrack,basemodfreqTrack,genesTrack,chrTrack])
        close();
    }

    function select_wgs_chr(){
        setGeneralParams((params)=> {return {...params,layout:"horizontal"}});
        setRegionsList([{id:crypto.randomUUID(),chr:"",start:"",end:"",color:"#f4a460"}]);
        let svTrack= {id:crypto.randomUUID(),type:"sv",...defaultTrackValues.sv};
        let copynumberTrack= {id:crypto.randomUUID(),type:"copynumber",...defaultTrackValues.copynumber,margin_above:0};
        let chrTrack={id:crypto.randomUUID(),type:"chr_axis",...defaultTrackValues["chr_axis"],margin_above:0};
        setHighlightsList([]);
        setTracksList([svTrack,copynumberTrack,chrTrack])
        close();
    }

    function select_wgs_circos(){
        setGeneralParams((params)=> {return {...params,layout:"circular"}});

        const regions=[];
        const colors=["#98671F","#65661B","#969833","#CE151D","#FF1A25","#FF0BC8","#FFCBCC","#FF9931","#FFCC3A","#FCFF44","#C4FF40","#00FF3B",
        "#2F7F1E","#2800C6","#6A96FA","#98CAFC","#00FEFD","#C9FFFE","#9D00C6","#D232FA","#956DB5","#5D5D5D","#989898","#CBCBCB"];
        for (let i = 1; i < 23; i++) {
            regions.push({id:crypto.randomUUID(),"chr":i.toString(),"start":"","end":"",color:colors[i-1]});
        }
        regions.push({id:crypto.randomUUID(),"chr":"X","start":"","end":"",color:colors[22]});
        regions.push({id:crypto.randomUUID(),"chr":"Y","start":"","end":"",color:colors[23]});
        setRegionsList(regions);

        let svTrack= {id:crypto.randomUUID(),type:"sv",...defaultTrackValues.sv};
        let copynumberTrack= {id:crypto.randomUUID(),type:"copynumber",...defaultTrackValues.copynumber,margin_above:0};
        let chrTrack={id:crypto.randomUUID(),type:"chr_axis",...defaultTrackValues["chr_axis"],margin_above:0};
        setHighlightsList([]);
        setTracksList([svTrack,copynumberTrack,chrTrack])
        close();
    }

    return(
    <div className="TemplatePanel" ref={popover}>
        <h2> Select template </h2>
        <div className="TrackTypeGrid">
            <TemplateButton name="bigwig" onClick={select_bigwig} />
            <TemplateButton name="hic" onClick={select_hic}/>
            <TemplateButton name="asm" onClick={select_asm}/>
            <TemplateButton name="wgs_chr" onClick={select_wgs_chr} />
            <TemplateButton name="wgs_circos" onClick={select_wgs_circos}/>
        </div>
        <button onClick={close}>Cancel</button>
    </ div>
    )
}

function TemplateButton({name,onClick}){
    return(
        <div className="TemplateButton" onClick={onClick}>
            {name}

        </div>
    )
}
