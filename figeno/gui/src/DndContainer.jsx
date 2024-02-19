import React from 'react';
import "./style.css";
import { DndContext, closestCenter , DragOverlay, MeasuringStrategy} from '@dnd-kit/core';
import { arrayMove, SortableContext, verticalListSortingStrategy } from '@dnd-kit/sortable';

const measuringConfig = {
    droppable: {
      strategy: MeasuringStrategy.Always,
    }
  };


export function DndContainer({children,header,items,setItems,show_active_item, className="TrackContainer"}) {
    const [activeId, setActiveId] = React.useState(null);


    function handleDragStart(event) {
        setActiveId(event.active.id);
    }
    function handleDragEnd(event){
        setActiveId(null);
        const {active,over} = event;
        if (active.id!==over.id){
            let activeIndex=0;
            let overIndex=0;
            for (let i=0; i<items.length;i++){
                if (items[i].id===active.id) {activeIndex=i;}
                if (items[i].id===over.id) {overIndex=i;}
            }
            setItems(arrayMove(items,activeIndex,overIndex));
        }
    }
    return (
        <div className={className}>
                {header}
            <div className="TrackContainerContent">
                <DndContext collisionDetection={closestCenter} onDragEnd={handleDragEnd} onDragStart={handleDragStart} measuring={measuringConfig}>
                    <SortableContext items={items} strategy={verticalListSortingStrategy} >
                        {children}
                    </SortableContext>
                    <DragOverlay zIndex={999}>
                        {activeId ? show_active_item(activeId) : null}
                    </DragOverlay>
                </DndContext>
            </div>

            
        </div >
    )
}