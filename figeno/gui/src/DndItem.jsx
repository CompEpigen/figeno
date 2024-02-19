import React from 'react';
import { useSortable } from '@dnd-kit/sortable';
import {CSS} from '@dnd-kit/utilities';
import "./style.css";



export function DndItem({children, id,copy_item,delete_item, className}) {
    const {attributes, listeners, setNodeRef, transform,transition, isDragging} = useSortable({id:id})
    let style = {transform: CSS.Transform.toString(transform), transition , opacity: isDragging ? 0.4 : undefined}

    
    return (
        <div className={className} ref={setNodeRef} style={style}>
           {/* Left side: handle... */}
            <div className="leftPanel">
                <div  {...attributes} {...listeners} className="trackHandle">
                <svg width="8mm" height="8mm" version="1.1" viewBox="0 0 20 20" xmlSpace="preserve" xmlns="http://www.w3.org/2000/svg"><g fill="#888888"><g strokeWidth="0"><rect x="6.5" y="9.3" width="7" height="1.4" ry=".99999" opacity=".97208"/><rect x="6.5" y="6.8" width="7" height="1.4" ry=".99999" opacity=".97208"/><rect x="6.5" y="11.8" width="7" height="1.4" ry=".99999" opacity=".97208"/></g><path d="m6.6198 4.6457c-0.14216 0.17344-0.14684 0.39141-0.013982 0.56836 0.13279 0.17695 0.38274 0.28711 0.65457 0.28711h5.499c0.27182 0 0.52178-0.11016 0.65456-0.28711 0.13279-0.17695 0.12653-0.39492-0.01398-0.56836l-2.7495-3.375c-0.1359-0.16757-0.37805-0.26953-0.64051-0.26953-0.26245 0-0.50303 0.10195-0.6405 0.26953z" strokeWidth=".01353"/><path d="m6.6102 15.354c-0.14216-0.17344-0.14684-0.39141-0.013982-0.56836 0.13279-0.17695 0.38274-0.28711 0.65457-0.28711h5.499c0.27182 0 0.52178 0.11016 0.65456 0.28711 0.13279 0.17695 0.12653 0.39492-0.01398 0.56836l-2.7495 3.375c-0.1359 0.16758-0.37805 0.26953-0.64051 0.26953-0.26245 0-0.50303-0.10195-0.6405-0.26953z" strokeWidth=".01353"/></g><ellipse cx="-3.66" cy="9.2286" rx=".0019184" ry=".18764" fill="#61718d" opacity=".97208" strokeWidth="0"/></svg>
                </div>
            </div>
            {/* Main part */}
            <div className="centerPanel">
                {children}
            </div>

            {/* Right part: delete... */}
            <div className="rightPanel">
                <button className="copyButton" onClick={copy_item}>
                <svg width="5mm" height="5.5mm" version="1.1" viewBox="0 0 5 5.5" xmlns="http://www.w3.org/2000/svg">
                    <rect width="5" height="5" opacity="0" strokeWidth=".20833"/>
                    <g transform="matrix(.25 0 0 .25 2.5 2.75)">
                    <rect x="-10" y="-11" width="20" height="20" opacity="0" strokeWidth=".83333"/>
                    <g transform="matrix(.41667 0 0 .41667 3.8456e-7 -1)">
                    <g strokeLinecap="round">
                        <g transform="translate(6.5,-2.5)">
                        <path transform="translate(-26.5,-17.5)" d="m14.5 33.5v-32h16.293l7.707 7.707v24.293z" fill="#f2faff"/>
                        <path transform="translate(-26.5,-17.5)" d="m30.586 2 7.414 7.414v23.586h-23v-31h15.586m0.414-1h-17v33h25v-25z" fill="#788b9c"/>
                        </g>
                        <g transform="translate(14.5,-14.5)">
                        <path transform="translate(-34.5,-5.5)" d="m30.5 9.5v-8h0.293l7.707 7.707v0.293z" fill="#fff"/>
                        <path transform="translate(-34.5,-5.5)" d="m31 2.414 6.586 6.586h-6.586v-6.586m0-1.414h-1v9h9v-1z" fill="#788b9c"/>
                        </g>
                        <g transform="translate(-6.5,2.5)">
                        <path transform="translate(-13.5,-22.5)" d="m1.5 38.5v-32h16.293l7.707 7.707v24.293z" fill="#fff"/>
                        <path transform="translate(-13.5,-22.5)" d="m17.586 7 7.414 7.414v23.586h-23v-31h15.586m0.414-1h-17v33h25v-25z" fill="#788b9c"/>
                        </g>
                        <g transform="translate(1.5,-9.5)">
                        <path transform="translate(-21.5,-10.5)" d="m17.5 14.5v-8h0.293l7.707 7.707v0.293z" fill="#fff"/>
                        <path transform="translate(-21.5,-10.5)" d="m18 7.414 6.586 6.586h-6.586v-6.586m0-1.414h-1v9h9v-1z" fill="#788b9c"/>
                        </g>
                    </g>
                    </g>
                    </g>
                </svg>
                </button>
                <button className="deleteButton" onClick={delete_item}>
                <svg width="4.5mm" height="4.9mm" version="1.1" viewBox="0 0 448 512" xmlns="http://www.w3.org/2000/svg">
                    <path d="M135.2 17.7L128 32H32C14.3 32 0 46.3 0 64S14.3 96 32 96H416c17.7 0 32-14.3 32-32s-14.3-32-32-32H320l-7.2-14.3C307.4 6.8 296.3 0 284.2 0H163.8c-12.1 0-23.2 6.8-28.6 17.7zM416 128H32L53.2 467c1.6 25.3 22.6 45 47.9 45H346.9c25.3 0 46.3-19.7 47.9-45L416 128z" fill="#e74c3c"/>
                    </svg>
                </button>
            </div>
        </div>
    )
}
