import {React, Fragment,useState,useRef, useCallback} from 'react';
import useClickOutside from "./useClickOutside";
import LoadingIcons from 'react-loading-icons'
import "./style.css";



export function LoadingScreen({errorMessage, setErrorMessage,setLoadingscreenActive}){
    function handleButton(){
        setErrorMessage("");
        setLoadingscreenActive(false);
    }

    return(
    <div className="LoadingScreen">
        Generating the figure...
        {(errorMessage==="")?(
        <LoadingIcons.SpinningCircles fill="#2980b9" speed="1.2" width="100px" height="100px" />
        ):""}

        {(errorMessage!=="")? (
            <>
            <div>An error has occured:</div>
            <p>{errorMessage}</p>
            <button onClick={handleButton}>Close</button>
            </>
        ):""}
    </ div>
    )
}
