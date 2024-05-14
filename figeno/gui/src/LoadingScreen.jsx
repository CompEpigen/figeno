import {React, Fragment,useState,useRef, useCallback} from 'react';
import {v4 as uuid4} from 'uuid';
import useClickOutside from "./useClickOutside";
import LoadingIcons from 'react-loading-icons'
import "./style.css";



export function LoadingScreen({message, setMessage,status,setStatus, warning, setWarning, setLoadingscreenActive}){
    function handleButton(){
        setStatus("");
        setMessage("");
        setWarning("");
        setLoadingscreenActive(false);
    }

    function openImage(){
        fetch("/open_image",{headers: {'Content-Type': 'application/json'}, body: JSON.stringify({file:message}),
        method:"POST"})
    }
    function openDir(){
        fetch("/open_dir",{headers: {'Content-Type': 'application/json'}, body: JSON.stringify({file:message}),
        method:"POST"})
    }

    const popover = useRef();
    useClickOutside(popover, ()=>{if (status!="") handleButton();});

    return(
    <div className="LoadingScreen" ref={popover}>
        {(status==="")?(
            <>
            <h2>Generating the figure...</h2>
            <LoadingIcons.SpinningCircles fill="#2980b9" speed="1.2" width="100px" height="100px" />
            </>
        ):""}

        {(status==="success")?(
            <>
            <h2>The figure was successfully generated. ✅</h2>
            <p style={{textAlign:"left",whiteSpace:"pre-wrap",backgroundColor:"white",fontSize:"1.2em"}}>
                The figure was saved to {message}.
            </p>

            {(warning!="")?(
                <>
                
                <div style={{overflowY: "auto", maxHeight:"400px", marginTop:"0px"}}>
                <p style={{textAlign:"left",whiteSpace:"pre-wrap",backgroundColor:"white", fontSize:"1.2em"}}>Warning ⚠️</p>
                    <p style={{textAlign:"left",whiteSpace:"pre-wrap",backgroundColor:"white"}}>{warning}</p>
                </div>
                </>
            ):""}


            <div style={{display:"flex", gap:"5px"}}>
            <button style={{fontSize:"1.3em"}} onClick={openImage}>Open image</button>
            <button style={{fontSize:"1.3em"}} onClick={openDir}>Open directory</button>
            <button style={{fontSize:"1.3em"}} onClick={handleButton}>Close</button>
            </div>
            </>
        ):""}

        {(status==="known_error")?(
            <>
            <h2>An error has occured. ❌</h2>
            <p style={{textAlign:"left",whiteSpace:"pre-wrap",backgroundColor:"white",fontSize:"1.2em"}}>{message}</p>
            <button style={{fontSize:"1.3em"}} onClick={handleButton}>Close</button>
            </>
        ):""}

        {(status==="unknown_error")?(
            <>
            <h2>An unknown error has occured. ❌</h2>
            <p style={{textAlign:"left",whiteSpace:"pre-wrap",backgroundColor:"white",fontSize:"1.2em"}}>
                The traceback below might help you solve the error. Otherwise, you can ask for help by <a href="https://github.com/CompEpigen/figeno/issues" target="_blank">creating an issue on the figeno GitHub repository</a>. Please provide a copy of the error message below, and of the config file that you used.
            </p>
            <div style={{overflowY: "auto", maxHeight:"400px"}}>
            <p style={{textAlign:"left",whiteSpace:"pre-wrap",backgroundColor:"white"}}>{message}</p>
            </div>
            <button style={{fontSize:"1.3em"}} onClick={handleButton}>Close</button>
            </>
        ):""}

    </ div>
    )
}




