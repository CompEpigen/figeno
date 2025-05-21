import React, {useState, useEffect} from 'react';
import {v4 as uuid4} from 'uuid';


import { GeneralContainer } from './GeneralContainer';
import { OutputContainer } from './OutputContainer';
import { RegionsContainer } from './RegionsContainer';
import { HighlightsContainer } from './HighlightsContainer';
import { TracksContainer, defaultTrackValues} from "./TracksContainer"
import { LoadingScreen } from './LoadingScreen';
import { GeneSelector } from './GeneSelector';
import { FileDialog } from './FileDialog';

import {ColorButton, ColorPanel} from "./ColorButton"
import { TrackTypePanel } from './TrackTypePanel';
import { TemplatePanel } from './TemplatePanel';
import "./style.css";

export default function App() {
  const [generalParams,setGeneralParams] = React.useState({"layout":"horizontal","reference":"hg19"});
  const [outputParams,setOutputParams] = React.useState({file:"",dpi:"400",width:"180"});
  const [regionsList,setRegionsList] = React.useState([{id:uuid4(),chr:"",start:"",end:"",color:"#f4a460"}]);
  const [highlightsList,setHighlightsList] = React.useState([]);
  const [tracksList,setTracksList] = React.useState([]);

  const [colorPanelActive,setColorPanelActive] = useState(false);
  const [color,setColor] = useState("#333333");
  const [updateColor, setUpdateColor] = useState(()=>{return false});

  const [trackTypePanelActive,setTrackTypePanelActive] = useState(false);
  const [setTrackType, setSetTrackType] = useState(()=>{return false});

  const [geneSelectorActive,setGeneSelectorActive] = useState(false);
  const [updateRegion, setUpdateRegion] = useState(()=>{return false});

  const [templatePanelActive,setTemplatePanelActive] = useState(false);

  const [loadingscreenActive,setLoadingscreenActive] = useState(false);
  const [message,setMessage] = useState("");
  const [warning,setWarning] = useState("");
  const [status,setStatus] = useState("");

  const [fileDialogActive,setFileDialogActive] = useState(false);
  const [fileDialogData,setFileDialogData] = useState({current_dir:"",dirs:[],files:[],file:"",dialog_type:"",update:function() {return}})

  let show_regions_color=false;
  for (const t of tracksList){
    if (t.type=="chr_axis" && (t.style=="arrow" || generalParams.layout=="circular")) show_regions_color=true; 
  }
  

  function openColorPanel(current_color,update_color){
    setColor(current_color);
    setUpdateColor(()=>update_color);
    setColorPanelActive(true);
  }

  function render_colorPanel(){
    if (colorPanelActive){
      return(
        <ColorPanel color={color} setColor={setColor} updateColor={updateColor} close={()=>{setColorPanelActive(false)}}></ColorPanel>
      )
    }
  }

  function openTrackTypePanel(set_type){
    setSetTrackType(()=>set_type);
    setTrackTypePanelActive(true);
  }

  function openGeneSelector(set_region){
    setUpdateRegion(()=>set_region);
    setGeneSelectorActive(true);
  }


  function set_regions_chromosomes(chromosomes){
    const regions=[];
      const colors=["#98671F","#65661B","#969833","#CE151D","#FF1A25","#FF0BC8","#FFCBCC","#FF9931","#FFCC3A","#FCFF44","#C4FF40","#00FF3B",
        "#2F7F1E","#2800C6","#6A96FA","#98CAFC","#00FEFD","#C9FFFE","#9D00C6","#D232FA","#956DB5","#5D5D5D","#989898","#CBCBCB"];
        for (let i = 0; i < chromosomes.length; i++) {
            regions.push({id:uuid4(),"chr":chromosomes[i].toString(),"start":"","end":"",color:colors[i%colors.length]});
        }
      setRegionsList(regions);
  }
  function add_all_chromosomes(){
    if (generalParams.reference==="hg19" || generalParams.reference==="hg38" || generalParams.reference==="mm10"){
      const n = (generalParams.reference==="mm10")?20:23;
      const chromosomes=[];
      for (let i = 1; i < n; i++) chromosomes.push(i.toString());
      chromosomes.push("X");
      chromosomes.push("Y");
      set_regions_chromosomes(chromosomes);
      
    }
    else{
      fetch("/get_all_chromosomes",{headers: {'Content-Type': 'application/json'}, body: JSON.stringify(generalParams),
        method:"POST"}).then(res => res.json()).then((data)=>{
          if (data.status!="success"){
            alert(data.message);
            const chromosomes=[];
            for (let i = 1; i < 23; i++) chromosomes.push(i.toString());
            chromosomes.push("X");
            chromosomes.push("Y");
            set_regions_chromosomes(chromosomes);
          }
          else set_regions_chromosomes(data.chromosomes);
        });
    }
    
  }

  function get_config(){
    const regions_out=[];
    for (const region of regionsList){
      const reg={chr:region.chr};
      if (region.hasOwnProperty("start") && (region.start!==null)){
        if (Number.isInteger(region.start)) {reg.start=region.start}
        else if (region.start.length>0) {reg.start=parseInt(region.start.replaceAll(",",""))}
      };
      if (region.hasOwnProperty("end") && (region.end!==null)){
        if (Number.isInteger(region.end)) {reg.end=region.end}
        else if (region.end.length>0) {reg.end=parseInt(region.end.replaceAll(",",""))}
      };
      if (show_regions_color){reg.color = region.color;}
      regions_out.push(reg);
    }

    const highlights_out=[];
    for (const region of highlightsList){
      const reg={chr:region.chr};
      if (region.hasOwnProperty("start") && (region.start!==null)){
        if (Number.isInteger(region.start)) {reg.start=region.start}
        else if (region.start.length>0) {reg.start=parseInt(region.start.replaceAll(",",""))}
      };
      if (region.hasOwnProperty("end") && (region.end!==null)){
        if (Number.isInteger(region.end)) {reg.end=region.end}
        else if (region.end.length>0) {reg.end=parseInt(region.end.replaceAll(",",""))}
      };
      reg.color=region.color;
      reg.opacity=parseFloat(region.opacity);
      highlights_out.push(reg);
    }

    const tracks_out=[];
    for (const track of tracksList){
      const t = {...track,height:parseFloat(track.height),margin_above:parseFloat(track.margin_above),fontscale:parseFloat(track.fontscale)};
      
      if (t.type=="chr_axis"){
        if (t.ticklabels_pos=="none"){
          delete t.unit;
          delete t.ticks_interval;
          delete t.ticks_angle;
        }
      }
      else if (t.type=="basemod_freq"){
        t.bams = t.bams.map((b)=>{const newBam={...b};delete newBam.id; return newBam;})
        t.bedmethyls = t.bedmethyls.map((b)=>{const newBed={...b};delete newBed.id; return newBed;})
      }
      else if (t.type=="bigwig" || t.type=="coverage"){
        if (t.scale!="custom"){
          delete t.scale_max;
          
        }
        if (t.scale!="group auto" && t.scale!="group auto per region"){
          delete t.group;
        }
        if (!t.show_negative && t.type=="bigwig"){
          delete t.negative_color;
        }
      }
      else if (t.type=="hic"){
        if (t.scale=="auto"){
          delete t.scale_min;
          delete t.scale_max;
        }
        else{
          delete t.scale_max_percentile;
        }
      }
      else if (t.type=="alignments"){
        if (t.group_by=="none"){
          delete t.show_unphased;
          delete t.exchange_haplotypes;
          delete t.show_haplotype_colors;
          delete t.haplotype_colors;
          delete t.haplotype_labels;
        }
        if (t.color_by=="none"){
          delete t.color_unmodified;
          delete t.basemods;
          delete t.fix_hardclip_basemod;
        }
        if (!t.link_splitreads){
          delete t.splitread_color;
          delete t.link_color;
          delete t.link_lw;
          delete t.only_show_splitreads;
          delete t.min_splitreads_breakpoints;
        }
        
      }
      else if (t.type=="copynumber"){
        if (t.input_type!="freec"){
          if (t.hasOwnProperty("freec_ratios")) delete t.freec_ratios;
          if (t.hasOwnProperty("freec_CNAs")) delete t.freec_CNAs;
        }
        if (t.input_type!="purple"){
          if (t.hasOwnProperty("purple_cn")) delete t.purple_cn;
          if (t.hasOwnProperty("color_cnloh")) delete t.color_cnloh;
        }
        if (t.input_type!="delly"){
          if (t.hasOwnProperty("delly_cn")) delete t.delly_ratios;
          if (t.hasOwnProperty("delly_CNAs")) delete t.delly_CNAs;
        }
      }
      delete t.id;
      tracks_out.push(t);
    }
    const config_dict = {
      general:generalParams,
      output: {...outputParams,dpi:parseInt(outputParams.dpi),width:parseInt(outputParams.width)},
      regions: regions_out,
      highlights: highlights_out,
      tracks: tracks_out
    };
    return config_dict;
  }

  function save_config(){
    const config_dict = get_config();

    function save_configfile(path){
      fetch("/save_config",{headers: {'Content-Type': 'application/json'}, body: JSON.stringify({path:path,config:config_dict}),
        method:"POST"})
    }
    

    fetch('/browse',{headers: {'Content-Type': 'application/json'}, body: JSON.stringify({path:"",dialog_type:"save_config"}),
    method:"POST"}).then(res => res.json()).then(data => {
      if (data.hasOwnProperty("current_dir")){
        setFileDialogData({current_dir:data.current_dir,dirs:data.dirs,files:data.files,file:"config.json",dialog_type:"save_config",update:save_configfile})
        setFileDialogActive(true);
      }
      else{
        save_configfile(data.path)
      }
      
      });
  }

  function run_figeno(){
    if (loadingscreenActive) return;
    const config_dict = get_config();
    setLoadingscreenActive(true);
    fetch("/run",{headers: {'Content-Type': 'application/json'}, body: JSON.stringify(config_dict),
        method:"POST"}).then(res => res.json()).then((data)=>{
          setStatus(data.status);
          setMessage(data.message);
          setLoadingscreenActive(true);
          if (data.hasOwnProperty("warning")) setWarning(data.warning);
        });
  }

  function load_config(){
    function load_data(data){
      if (!data.hasOwnProperty("general")){return}
      setGeneralParams(data.general);
      setOutputParams(data.output);
      
      const regions=[];
      for (const r of data.regions){
        const region={id:uuid4(),...r};
        if (!region.hasOwnProperty("color")) {region.color="#f4a460"}
        regions.push(region);
      }
      setRegionsList(regions);

      if (data.hasOwnProperty("highlights")){
        const highlights=[];
        for (const r of data.highlights){
          const region={id:uuid4(),...r};
          highlights.push(region);
        }
      setHighlightsList(highlights);
      }

      const tracks=[];
      for (const t of data.tracks){
        if (t.hasOwnProperty("type")){
          const track={id:uuid4(),type:t.type,...defaultTrackValues[t.type],...t};

          if (t.type=="basemod_freq"){
            let bams = [];
            if (t.hasOwnProperty("bams")){
              for (const x of t.bams){
                bams.push({id:uuid4(),...x})
              }
            }
            track.bams = bams;
            let bedmethyls = [];
            if (t.hasOwnProperty("bedmethyls")){
              for (const x of t.bedmethyls){
                bedmethyls.push({id:uuid4(),...x})
              }
            }
            track.bedmethyls = bedmethyls;
          }
          tracks.push(track);
        }
      }
      setTracksList(tracks);
    }

    function load_configfile(path){
      fetch('/load_config',{headers: {'Content-Type': 'application/json'}, body: JSON.stringify({path:path}),
      method:"POST"}).then(res => res.json()).then(data => {
        load_data(data)
      })
    }


    
    fetch('/browse',{headers: {'Content-Type': 'application/json'}, body: JSON.stringify({path:fileDialogData.current_dir,dialog_type:"load_config"}),
    method:"POST"}).then(res => res.json()).then(data => {
      if (data.hasOwnProperty("current_dir")){
        setFileDialogData({current_dir:data.current_dir,dirs:data.dirs,files:data.files,file:"",dialog_type:"load_config",update:load_configfile})
        setFileDialogActive(true);
      }
      else{
        load_configfile(data.path)
      }
      
      });

    /*fetch('/load_config').then(res => res.json()).then(data => {
      
    });*/
  }

  return (
    <div className="app">
    
    <div className='topBar'>
      <button onClick={save_config} className="buttonImage blueButton">
        <svg xmlns="http://www.w3.org/2000/svg" height="24" viewBox="0 -960 960 960" width="24"><path d="M800-663.077v438.462Q800-197 781.5-178.5 763-160 735.385-160h-510.77Q197-160 178.5-178.5 160-197 160-224.615v-510.77Q160-763 178.5-781.5 197-800 224.615-800h438.462L800-663.077ZM760-646 646-760H224.615q-10.769 0-17.692 6.923T200-735.385v510.77q0 10.769 6.923 17.692T224.615-200h510.77q10.769 0 17.692-6.923T760-224.615V-646ZM480-298.461q33.077 0 56.539-23.462Q560-345.384 560-378.461T536.539-435Q513.077-458.462 480-458.462T423.461-435Q400-411.538 400-378.461t23.461 56.538q23.462 23.462 56.539 23.462Zm-209.231-270.77h296.924v-120H270.769v120ZM200-646v446-560 114Z" fill="#ffffff"/></svg>
          Save config
      </button>
          
        <button onClick={load_config} className="buttonImage blueButton">
        <svg width="4.5mm" height="4.4mm" version="1.1" viewBox="0 0 576 512" xmlns="http://www.w3.org/2000/svg">
                  <path d="M384 480h48c11.4 0 21.9-6 27.6-15.9l112-192c5.8-9.9 5.8-22.1 .1-32.1S555.5 224 544 224H144c-11.4 0-21.9 6-27.6 15.9L48 357.1V96c0-8.8 7.2-16 16-16H181.5c4.2 0 8.3 1.7 11.3 4.7l26.5 26.5c21 21 49.5 32.8 79.2 32.8H416c8.8 0 16 7.2 16 16v32h48V160c0-35.3-28.7-64-64-64H298.5c-17 0-33.3-6.7-45.3-18.7L226.7 50.7c-12-12-28.3-18.7-45.3-18.7H64C28.7 32 0 60.7 0 96V416c0 35.3 28.7 64 64 64H87.7 384z" fill="#ffffff"/>
          </svg>
          Load config
      </button>

      <button onClick={()=>setTemplatePanelActive(true)} className="buttonImage blueButton">
      <svg xmlns="http://www.w3.org/2000/svg" height="25" viewBox="0 -960 960 960" width="25"><path d="M340-240h435.385q9.23 0 16.923-7.692Q800-255.385 800-264.615V-367H340v127ZM160-593h140v-127H184.615q-9.23 0-16.923 7.692Q160-704.615 160-695.385V-593Zm0 187h140v-147H160v147Zm24.615 166H300v-127H160v102.385q0 9.23 7.692 16.923Q175.385-240 184.615-240ZM340-406h460v-147H340v147Zm0-187h460v-102.385q0-9.23-7.692-16.923Q784.615-720 775.385-720H340v127ZM184.615-200Q157-200 138.5-218.5 120-237 120-264.615v-430.77Q120-723 138.5-741.5 157-760 184.615-760h590.77Q803-760 821.5-741.5 840-723 840-695.385v430.77Q840-237 821.5-218.5 803-200 775.385-200h-590.77Z" fill="#ffffff"/></svg>
          Open template
      </button>

      <button onClick={run_figeno} className="buttonImage greenButton">
      <svg width="12" height="14.5" version="1.1" viewBox="0 0 384 464" xmlns="http://www.w3.org/2000/svg">
        <path d="m73 13.458c-14.8-9.1-33.4-9.4-48.5-0.9-15.1 8.5-24.5 24.5-24.5 41.9v352c0 17.4 9.4 33.4 24.5 41.9s33.7 8.1 48.5-0.9l288-176c14.3-8.7 23-24.2 23-41s-8.7-32.2-23-41z" fill="#ffffff"/>
      </svg>
        Generate figure
        </button>
    </div>
    {render_colorPanel()}
    {fileDialogActive ? (<FileDialog fileDialogData={fileDialogData} setFileDialogData={setFileDialogData} close={()=>setFileDialogActive(false)}/>):""}
    {trackTypePanelActive ? (<TrackTypePanel setTrackType={setTrackType} close={()=>setTrackTypePanelActive(false)}/>):""}
    {geneSelectorActive ? (<GeneSelector updateRegion={updateRegion} setGeneSelectorActive={setGeneSelectorActive} generalParams={generalParams}/>):""}
    {loadingscreenActive ? (<LoadingScreen message={message} setMessage={setMessage} status={status} setStatus={setStatus} warning={warning} setWarning={setWarning} setLoadingscreenActive={setLoadingscreenActive} />):""}
    {templatePanelActive && (<TemplatePanel setTracksList={setTracksList} setRegionsList={setRegionsList} setHighlightsList={setHighlightsList} setGeneralParams={setGeneralParams} close={()=>setTemplatePanelActive(false)}/>)}
      <div style={{display:"grid",gridTemplateColumns:"repeat(auto-fit,minmax(520px,1fr))", width:"100%", gap:"10px"}}>
      <GeneralContainer generalParams={generalParams} setGeneralParams={setGeneralParams}/>
      <OutputContainer outputParams={outputParams} setOutputParams={setOutputParams} setFileDialogData={setFileDialogData} setFileDialogActive={setFileDialogActive}/>
      </div>

      <div style={{display:"grid",gridTemplateColumns:"repeat(auto-fit,minmax(520px,1fr))", width:"100%", gap:"10px"}}>
      <RegionsContainer key="regions" regionsList={regionsList} setRegionsList={setRegionsList} openColorPanel={openColorPanel} openGeneSelector={openGeneSelector} add_all_chromosomes={add_all_chromosomes} show_regions_color={show_regions_color}/>
      <HighlightsContainer key="highlights" regionsList={highlightsList} setRegionsList={setHighlightsList} openColorPanel={openColorPanel}/>
      </div>


        <div>
          <TracksContainer key="tracks" tracksList={tracksList} setTracksList={setTracksList} openColorPanel={openColorPanel} openTrackTypePanel={openTrackTypePanel} setFileDialogData={setFileDialogData} setFileDialogActive={setFileDialogActive} />
        </div>

    </div>
   
  );
  
}