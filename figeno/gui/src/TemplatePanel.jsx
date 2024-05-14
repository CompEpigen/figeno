import {React, useRef} from 'react';
import {v4 as uuid4} from 'uuid';
import useClickOutside from "./useClickOutside";
import "./style.css";
import { defaultTrackValues } from './TracksContainer';

import bigwig_figure from './images/template_bigwig.png';
import hic_figure from './images/template_hic.png';
import asm_figure from './images/template_asm.png';
import wgs_figure from './images/template_wgs_chr.png';
import circos_figure from './images/template_wgs_circos.png';


export function TemplatePanel({setTracksList,setRegionsList,setHighlightsList,setGeneralParams,close}){
    const templates = ["bigwig", "hic", "asm", "wgs_chr", "wgs_circos"];
    const popover = useRef();
    useClickOutside(popover, close);

    function select_bigwig(){
        setGeneralParams((params)=> {return {...params,layout:"horizontal"}});
        setRegionsList([{id:uuid4(),chr:"",start:"",end:"",color:"#f4a460"}]);
        let bigwigTrack= {id:uuid4(),type:"bigwig",...defaultTrackValues.bigwig};
        let genesTrack={id:uuid4(),type:"genes",...defaultTrackValues.genes};
        let chrTrack={id:uuid4(),type:"chr_axis",...defaultTrackValues["chr_axis"]};
        setHighlightsList([]);
        setTracksList([bigwigTrack,genesTrack,chrTrack])
        close();
    }

    function select_hic(){
        setGeneralParams((params)=> {return {...params,layout:"horizontal"}});
        setRegionsList([{id:uuid4(),chr:"",start:"",end:"",color:"#f4a460"}]);
        let hicTrack= {id:uuid4(),type:"hic",...defaultTrackValues.hic};
        let genesTrack={id:uuid4(),type:"genes",...defaultTrackValues.genes};
        let chrTrack={id:uuid4(),type:"chr_axis",...defaultTrackValues["chr_axis"]};
        setHighlightsList([]);
        setTracksList([hicTrack,genesTrack,chrTrack])
        close();
    }

    function select_asm(){
        setGeneralParams((params)=> {return {...params,layout:"horizontal"}});
        setRegionsList([{id:uuid4(),chr:"",start:"",end:"",color:"#f4a460"}]);
        let alignmentsTrack= {id:uuid4(),type:"alignments",...defaultTrackValues.alignments,group_by:"haplotype",color_by:"basemod"};
        let basemodfreqTrack= {id:uuid4(),type:"basemod_freq",...defaultTrackValues["basemod_freq"],
                        bams:[{id:uuid4(),"file": "","base": "C","mod": "m","min_coverage": 6,"linewidth": 3,"opacity": 1,
                        "fix_hardclip": false,"split_by_haplotype": true,"colors": ["#27ae60","#e67e22"],"labels":["",""]}]};
        let genesTrack={id:uuid4(),type:"genes",...defaultTrackValues.genes};
        let chrTrack={id:uuid4(),type:"chr_axis",...defaultTrackValues["chr_axis"]};
        setHighlightsList([]);
        setTracksList([alignmentsTrack,basemodfreqTrack,genesTrack,chrTrack])
        close();
    }

    function select_wgs_chr(){
        setGeneralParams((params)=> {return {...params,layout:"horizontal"}});
        setRegionsList([{id:uuid4(),chr:"",start:"",end:"",color:"#f4a460"}]);
        let svTrack= {id:uuid4(),type:"sv",...defaultTrackValues.sv};
        let copynumberTrack= {id:uuid4(),type:"copynumber",...defaultTrackValues.copynumber,margin_above:0};
        let chrTrack={id:uuid4(),type:"chr_axis",...defaultTrackValues["chr_axis"],margin_above:0,unit:"Mb"};
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
            regions.push({id:uuid4(),"chr":i.toString(),"start":"","end":"",color:colors[i-1]});
        }
        regions.push({id:uuid4(),"chr":"X","start":"","end":"",color:colors[22]});
        regions.push({id:uuid4(),"chr":"Y","start":"","end":"",color:colors[23]});
        setRegionsList(regions);

        let svTrack= {id:uuid4(),type:"sv",...defaultTrackValues.sv};
        let copynumberTrack= {id:uuid4(),type:"copynumber",...defaultTrackValues.copynumber,margin_above:0,max_cn:3.9,grid_major:false,grid_minor:false};
        let chrTrack={id:uuid4(),type:"chr_axis",...defaultTrackValues["chr_axis"],margin_above:0};
        setHighlightsList([]);
        setTracksList([svTrack,copynumberTrack,chrTrack])
        close();
    }

    return(
    <div className="TemplatePanel" ref={popover}>
        <h2> Select template </h2>
        <div className="TrackTypeGrid">
            <TemplateButton name="bigwig" onClick={select_bigwig}><img src={bigwig_figure} width="200" height="40" alt="hic" /></TemplateButton>
            <TemplateButton name="hic" onClick={select_hic}><img src={hic_figure} width="200" height="100" alt="hic" /></TemplateButton>
            <TemplateButton name="asm" onClick={select_asm}><img src={asm_figure} width="200" height="100" alt="hic" /></TemplateButton>
            <TemplateButton name="wgs_chr" onClick={select_wgs_chr} ><img src={wgs_figure} width="200" height="50" alt="hic" /></TemplateButton>
            <TemplateButton name="wgs_circos" onClick={select_wgs_circos}><img src={circos_figure} width="100" height="100" alt="hic" /></TemplateButton>
        </div>
        <button onClick={close} style={{fontSize:"1.2em", marginTop:"8px"}}>Cancel</button>
    </ div>
    )
}

function TemplateButton({children,name,onClick}){
    return(
        <div className="TemplateButton" onClick={onClick}>
            {name}
            {children}

        </div>
    )
}
