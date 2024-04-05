import React from 'react';
import "./style.css";
import { PathEntry } from './Track';


export function GeneralContainer({generalParams,setGeneralParams}) {
    
    function set_value(attribute,value){
        const newVal={...generalParams};
        newVal[attribute]= value;
        setGeneralParams(newVal);
    }

    function handle_reference_change(value){
        if (value!=generalParams.reference){
            const newVal={...generalParams,reference:value};
            if (value=="custom"){
                newVal["genes_file"] = "";
                newVal["chrarms_file"] = "";
                newVal["cytobands_file"] = "";

            }
            else if (generalParams.reference=="custom"){
                delete newVal.genes_file;
                delete newVal.chrarms_file;
                delete newVal.cytobands_file;
            }
            setGeneralParams(newVal);
        }
    }
    




    const header=<div className="TrackContainerHeader">
        <h3>General</h3>
    </div>;

    return (
        <div className="TrackContainer">
            {header}
            <div className="TrackContainerContent">
                <div className="trackOption" style={{paddingLeft:"15px"}}>
                    <div className='formItem'>
                        <label htmlFor={"layout"}>Layout:</label>
                        <select id={"layout"} value={generalParams.layout} onChange={(e) =>{set_value("layout",e.target.value)}}> 
                                    <option className="dropDownOption" key="horizontal" value="horizontal">horizontal</option>
                                    <option className="dropDownOption" key="circular" value="circular">circular</option>
                                    <option className="dropDownOption" key="symmetrical" value="symmetrical">symmetrical</option>
                                    <option className="dropDownOption" key="stacked" value="stacked">stacked</option>
                        </select>
                    </div>
                    <div className='formItem'>
                        <label htmlFor={"reference"}>Reference:</label>
                        <select id={"reference"} value={generalParams.reference} onChange={(e) =>{handle_reference_change(e.target.value)}}> 
                                    <option className="dropDownOption" key="hg19" value="hg19">hg19</option>
                                    <option className="dropDownOption" key="hg38" value="hg38">hg38</option>
                                    <option className="dropDownOption" key="custom" value="custom">custom</option>
                        </select>
                    </div>

                    {generalParams.reference==="custom"?(<>
                        <PathEntry id={"genes_file"} label="Genes file:" value={generalParams.genes_file} set_value={(val) => set_value("genes_file",val)}/>
                        <PathEntry id={"chrarms_file"} label="Chr arms file:" value={generalParams.chrarms_file} set_value={(val) => set_value("chrarms_file",val)}/>
                        <PathEntry id={"cytobands"} label="Cytobands file:" value={generalParams.cytobands_file} set_value={(val) => set_value("cytobands_file",val)}/>
                    </>):""}
                </div>
            </div>
        </div>
    )
}