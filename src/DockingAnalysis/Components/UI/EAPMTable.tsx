import { AgGridReact } from "ag-grid-react"; // the AG Grid React Component

import { useCallback, useRef } from "react";
import { getFileForMolstar, postData } from "../../Common/utils";

function parseColumns(columns: Array<string>): Array<any> {
  const parsed: Array<any> = [];
  columns.forEach((c) => {
    parsed.push({ field: c });
  });

  return parsed;
}

export function EAPMTable({ results }: { results: any }) {
  const gridRef = useRef<AgGridReact>();

  const cellClickedListener = useCallback(async (event: any) => {
    // If the user clicked on the show column
    if (/^(ligand|system|complex)$/.test(event.colDef.field.toLowerCase())) {
      //   const postTO = window.location.href + "api/getPDB";
      //   const data = await postData(postTO, {
      //     path: block.path,
      //     directPDB: event.value,
      //   });
      //   await parent.molstar.loadMoleculeFile(
      //     getFileForMolstar({ contents: data.data, filename: data.filename })
      //   );
    }
  }, []);

  return (
    <>
      <div className="ag-theme-alpine w-full">
        <AgGridReact
          // @ts-ignore
          className="text-center justify-center items-center"
          // @ts-ignore
          ref={gridRef} // Ref for accessing Grid's API
          rowData={results.rows} // Row Data for Rows
          columnDefs={parseColumns(results.columns)} // Column Defs for Columns
          defaultColDef={{
            flex: 1,
            resizable: true,
            sortable: true,
            filter: true,
            floatingFilter: true,
            floatingFilterComponentParams: {
              suppressFilterButton: true,
            },
            minWidth: 150,
            // Aing the text to the center
            cellStyle: { textAlign: "center" },
          }} // Default Column Properties
          animateRows={true} // Optional - set to 'true' to have rows animate when sorted
          //   onCellClicked={cellClickedListener} // Optional - registering for Grid Event
          getRowId={(params) => params.data.id}
          domLayout="autoHeight"
        />
      </div>
    </>
  );
}
