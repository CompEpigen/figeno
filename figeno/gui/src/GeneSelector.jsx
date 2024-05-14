import {React, Fragment,useState, useRef} from 'react';
import useClickOutside from './useClickOutside';
import "./style.css";

export function GeneSelector({ updateRegion,setGeneSelectorActive, generalParams}){
    const popover = useRef();
    useClickOutside(popover, ()=>setGeneSelectorActive(false));
    const [errorMessage,setErrorMessage]=useState("");
    const [status,setStatus]=useState("");
    const [genename,setGenename]=useState("");
    function handleButton(){
      setStatus("")
        setErrorMessage("");
        setGeneSelectorActive(false);
    }

    function handleClickOK(){
        let genes_file="";
        if (generalParams.hasOwnProperty("genes_file")) genes_file=generalParams.genes_file;
        fetch("/find_gene",{headers: {'Content-Type': 'application/json'}, body: JSON.stringify({gene_name:genename,reference:generalParams.reference,genes_file:genes_file}),
        method:"POST"}).then(res => res.json()).then((data)=>{
          if (data.status=="success"){
            updateRegion(data.chr,data.start,data.end);
            setGeneSelectorActive(false);
          }
          else{
            setStatus(data.status)
            setErrorMessage(data.message);
          }
        });
    }

    return(
    <div className="LoadingScreen GeneSelector" ref={popover}>
        <h3 style={{marginBottom:"15px"}}>Center the region around a gene.</h3>
        <div className='formItem'>
                <label htmlFor={"gene_name_selection"} style={{fontSize:"1.2em"}}>Gene name:</label>
                <input autoFocus id={"gene_name_selection"} style={{width:"8em", fontSize:"1.3em",margin:"5px"}} value={genename} onKeyDown={(e)=>{if (e.key=="Enter") handleClickOK();}} onChange={(e) => setGenename(e.target.value)}/>
                <button onClick={handleClickOK} style={{fontSize:"1.2em"}}>OK</button>
        </div>
        
        

        {(status=="known_error")? (
            <>
            <div style={{fontSize:"1.2em", marginTop:"10px"}} >An error has occured:</div>
            <div style={{overflowY: "auto", maxHeight:"400px"}}>
            <p style={{textAlign:"left",whiteSpace:"pre-wrap",backgroundColor:"white"}}>{errorMessage}</p>
            </div>
            </>
        ):""}
        {(status=="unknown_error")? (
            <>
            <div style={{fontSize:"1.2em", marginTop:"10px"}}>An error has occured:</div>
            <div style={{overflowY: "auto", maxHeight:"400px"}}>
            <p style={{textAlign:"left",whiteSpace:"pre-wrap",backgroundColor:"white"}}>{errorMessage}</p>
            </div>
            </>
        ):""}
    <button onClick={handleButton} style={{margin:"10px",fontSize:"1.2em"}}>Close</button>
    </ div>
    )
}