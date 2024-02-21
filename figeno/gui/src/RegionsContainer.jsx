import React from 'react';
import {v4 as uuid4} from 'uuid';
import "./style.css";
import {Region} from './Region'
import { DndContainer } from './DndContainer';


export function RegionsContainer({regionsList,setRegionsList,openColorPanel,show_regions_color}) {
    

    function set_value(id,attribute,value){
        setRegionsList((regions)=>{
            return regions.map((region)=>{
                if (region.id==id) {
                    region[attribute]=value;
                    return region;}
                else {return region;}
            })
        });
    }
    


    function add_region(){
        setRegionsList([...regionsList,{id:uuid4(),chr:"",start:"",end:"",color:"#f4a460"}]);
    }
    function delete_region(id){
        setRegionsList((regions_list) => {return regions_list.filter(region => region.id!==id)})
    }
    function copy_region(id){
        const index = regionsList.findIndex((t)=>t.id==id);
        const newTrack = {...regionsList[index],id:uuid4()};
        const newList = [...regionsList]
        newList.splice(index,0,newTrack);
        setRegionsList(newList);
    }
    function getRegionById (regions, id) {
        return regions.find((region) => region.id === id);
    }
    function add_all_chromosomes(){
        const regions=[];
        const colors=["#98671F","#65661B","#969833","#CE151D","#FF1A25","#FF0BC8","#FFCBCC","#FF9931","#FFCC3A","#FCFF44","#C4FF40","#00FF3B",
        "#2F7F1E","#2800C6","#6A96FA","#98CAFC","#00FEFD","#C9FFFE","#9D00C6","#D232FA","#956DB5","#5D5D5D","#989898","#CBCBCB"];
        for (let i = 1; i < 23; i++) {
            regions.push({id:uuid4(),"chr":i.toString(),"start":"","end":"",color:colors[i-1]});
        }
        regions.push({id:uuid4(),"chr":"X","start":"","end":"",color:colors[22]});
        regions.push({id:uuid4(),"chr":"Y","start":"","end":"",color:colors[23]});
        setRegionsList(regions);


    }


    const header=<div className="TrackContainerHeader">
        <h3>Regions</h3>
        <button className="ButtonHeader" onClick={add_region}>
        <svg xmlns="http://www.w3.org/2000/svg" height="24" viewBox="0 -960 960 960" width="24"><path d="M440-440H200v-80h240v-240h80v240h240v80H520v240h-80v-240Z" fill="#ffffff" /></svg>
            Add region
        </button>
        <button className="ButtonHeader" onClick={add_all_chromosomes}>
        <svg xmlns="http://www.w3.org/2000/svg" height="24" viewBox="0 -960 960 960" width="24"><path d="M440-440H200v-80h240v-240h80v240h240v80H520v240h-80v-240Z" fill="#ffffff" /></svg>
            Add all chromosomes
        </button>
    </div>;

    function show_active_item(active_id){
        return <Region region={{...getRegionById(regionsList,active_id),id:active_id+"overlay"}} className={"track trackOverlay"} show_region_color={show_regions_color}/> 
    }

    return (
        <>
        <DndContainer header={header} items={regionsList} setItems={setRegionsList} show_active_item={show_active_item} >
        {regionsList.map((region)=>{
                    return <Region key={region.id} region={region} set_value={(attribute,value)=>set_value(region.id,attribute,value)}  className={"track"} copy_region={copy_region} delete_region={delete_region} openColorPanel={openColorPanel} show_region_color={show_regions_color}/>
                })}
        </DndContainer>
        </>
    )
}