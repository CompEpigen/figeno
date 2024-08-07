import React from 'react';
import {v4 as uuid4} from 'uuid';
import "./style.css";
import {Region} from './Region'
import { DndContainer } from './DndContainer';


export function RegionsContainer({regionsList,setRegionsList,openColorPanel, openGeneSelector,add_all_chromosomes,show_regions_color}) {
    
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
                    return <Region key={region.id} region={region} set_value={(attribute,value)=>set_value(region.id,attribute,value)}  className={"track"} copy_region={copy_region} delete_region={delete_region} openColorPanel={openColorPanel} openGeneSelector={openGeneSelector} show_region_color={show_regions_color}/>
                })}
        </DndContainer>
        </>
    )
}