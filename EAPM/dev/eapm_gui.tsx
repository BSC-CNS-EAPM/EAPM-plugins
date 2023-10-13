import React, {
  useMemo,
  useCallback,
  useEffect,
  useRef,
  useState,
} from "react";
import { horusGet } from "../../../../Views/Utils/utils";
import Plot, { PlotParams } from "react-plotly.js";
import NBDSuiteData from "./nbdsuitedata";
import HorusMolstar from "../../../../Views/Components/Molstar/HorusWrapper/horusmolstar";

// CSS style
import "./nbdsuite.css";

import { AgGridReact } from "ag-grid-react"; // the AG Grid React Component
import "ag-grid-community/styles/ag-grid.css"; // Core grid CSS, always needed
import "ag-grid-community/styles/ag-theme-alpine.css"; // Optional theme CSS

import { HorusModal } from "../../../../Views/Components/reusable";

declare global {
  interface Window {
    molstar?: HorusMolstar;
  }
}

interface OpenFolderProps {
  setOpenedFolder: (value: React.SetStateAction<string>) => void;
}

function OpenFolder(props: OpenFolderProps) {
  const openPickFolder = async () => {
    const request = await horusGet("/openfolder");

    const data = await request.json();

    const folder = data.path;

    props.setOpenedFolder(folder);
  };

  return (
    <div>
      <button onClick={openPickFolder}>Open Folder</button>
    </div>
  );
}

interface OpenInputFileProps {
  setOpenedFile: (value: React.SetStateAction<string>) => void;
}

function OpenInputFile(props: OpenInputFileProps) {
  const openPickFile = async () => {
    const request = await horusGet("/openfile");

    const data = await request.json();

    const file = data.path;

    props.setOpenedFile(file);
  };

  const openFolderIcon = (
    <div
      style={{
        height: "4rem",
        width: "4rem",
        marginBottom: "1rem",
      }}
    >
      <svg
        xmlns="http://www.w3.org/2000/svg"
        fill="none"
        viewBox="0 0 24 24"
        strokeWidth={1.5}
        stroke="currentColor"
        className="w-6 h-6"
      >
        <path
          strokeLinecap="round"
          strokeLinejoin="round"
          d="M3.75 9.776c.112-.017.227-.026.344-.026h15.812c.117 0 .232.009.344.026m-16.5 0a2.25 2.25 0 00-1.883 2.542l.857 6a2.25 2.25 0 002.227 1.932H19.05a2.25 2.25 0 002.227-1.932l.857-6a2.25 2.25 0 00-1.883-2.542m-16.5 0V6A2.25 2.25 0 016 3.75h3.879a1.5 1.5 0 011.06.44l2.122 2.12a1.5 1.5 0 001.06.44H18A2.25 2.25 0 0120.25 9v.776"
        />
      </svg>
    </div>
  );

  return (
    <div
      className="flex-center-col"
      style={{
        margin: "auto",
        justifyContent: "center",
        marginTop: "1rem",
        alignContent: "center",
        alignItems: "center",
        height: "100%",
      }}
    >
      {openFolderIcon}
      <button onClick={openPickFile}>Open NBDSuite input</button>
    </div>
  );
}

interface PELEPlotProps {
  nbdData: NBDSuiteData;
  addExtraData: (data: any) => void;
  updateSelectedComplex: (e: any) => void;
}

function PELEPlot(props: PELEPlotProps) {
  const nbdData = props.nbdData;

  const [xAxis, setXAxis] = useState("");
  const [yAxis, setYAxis] = useState("");

  const [data, setData] = useState<Array<any>>([]);
  const [legendData, setLegendData] = useState<Array<any>>([]);

  const [axisOptions, setAxisOptions] = useState<Array<string>>([]);

  useEffect(() => {
    const options = nbdData.axisOptions();

    setAxisOptions(options);

    setXAxis(options[0]);
    setYAxis(options[1]);
  }, []);

  function dataToPlot() {
    setData(nbdData.getPlot(xAxis, yAxis));

    // setLegendData(legendData);
  }

  useEffect(() => {
    dataToPlot();
  }, [xAxis, yAxis]);

  const updateXAxis = (e: any) => {
    setXAxis(e.target.value);
  };

  const updateYAxis = (e: any) => {
    setYAxis(e.target.value);
  };

  const xAxisOptions = (
    <select onChange={updateXAxis} key={"x-axis"}>
      {axisOptions.map((c: string) => {
        return <option key={c}>{c}</option>;
      })}
    </select>
  );

  const yAxisOptions = (
    <select onChange={updateYAxis} key={"y-axis"}>
      {axisOptions.map((c: string) => {
        return <option key={c}>{c}</option>;
      })}
    </select>
  );

  const complexesView = (
    <select key={"complex-select"}>
      {nbdData.complexes.map((c: any) => {
        return <option onChange={props.updateSelectedComplex}>{c}</option>;
      })}
    </select>
  );

  const ref = useRef();

  return (
    <div>
      <div
        className="blurred-squared-div"
        style={{
          overflow: "hidden",
        }}
      >
        <div>
          <div
            style={{
              display: "flex",
              justifyContent: "space-between",
              alignItems: "center",
            }}
          >
            <div className="title-color">{nbdData.simulationName}</div>
            <div
              style={{
                display: "flex",
                alignItems: "center",
                gap: "1rem",
              }}
            >
              <div>Select a complex:</div>
              {complexesView}
            </div>
          </div>
          <Plot
            ref={ref}
            key={"plot"}
            layout={{
              autosize: true,
              margin: {
                l: 50,
                r: 50,
                b: 50,
                t: 25,
                pad: 0,
              },
            }}
            data={data}
            onClick={(e) => {
              props.addExtraData(e.points[0].data.data);
            }}
            style={{
              width: "100%",
              height: "100%",
              margin: "0",
              padding: "0",
              overflow: "hidden",
            }}
          />
        </div>
        <div
          className="flex-center-row"
          style={{
            alignContent: "center",
            justifyContent: "center",
          }}
        >
          <div className="flex-center-row">
            {xAxisOptions}
            vs.
            {yAxisOptions}
          </div>
        </div>
      </div>
    </div>
  );
}

interface PELETableProps {
  nbdData: NBDSuiteData;
  extraData: Array<any>;
}

function PELETable(props: PELETableProps) {
  const { nbdData, extraData } = props;

  const gridRef = useRef<AgGridReact>(); // Optional - for accessing Grid's API
  const [rowData, setRowData] = useState<Array<any>>(); // Set rowData to Array of Objects, one Object per Row

  // Each Column Definition results in one Column.
  const [columnDefs, setColumnDefs] = useState<any>([]);

  // DefaultColDef sets props common to all Columns
  const defaultColDef = useMemo(() => ({
    sortable: true,
  }));

  // Example of consuming Grid Event
  const cellClickedListener = useCallback((event) => {
    // If the user clicked on the delete column
    if (event.colDef.headerName === "Delete") {
      const row = event.data;
      gridRef.current?.api.applyTransaction({ remove: [row] });
    }

    // If the user clicked on the show column
    if (event.colDef.headerName === "Add to view") {
      const row = event.data;
      nbdData.addPDBToMolstarFromRow(row);
    }

    // If the user clicked on the save column
    if (event.colDef.headerName === "Save") {
      const row = event.data;
      nbdData.savePDBfromRow(row);
    }
  }, []);

  function loadDataTable() {
    const { columns, rows } = nbdData.getTable();
    setColumnDefs(columns);
    setRowData(rows);
  }

  useEffect(() => {
    if (nbdData) {
      loadDataTable();
    }
  }, []);

  // Merge the extra data into the rowData
  useEffect(() => {
    if (extraData && rowData) {
      // Check if the row already exists in rowData
      const rowExists = gridRef.current?.api.getRowNode(extraData["id"]);

      if (!rowExists) {
        // Add to the extra data the delete column
        extraData["delete"] = "Delete";
        // Add star * to the data if its not coming from the repr
        if (extraData["Type"] !== "Representative") {
          extraData["Cluster label"] = extraData["Cluster label"] + "*";
        }
        gridRef.current?.api.applyTransaction({
          add: [extraData],
        });
        // setRowData([...rowData, extraData]);
      }
    }
  }, [extraData]);

  const [metricLoading, setMetricLoading] = useState(false);

  const newAtomAtomDistance = async () => {
    setMetricLoading(true);
    await nbdData.addAtomAtomDistance();
    setMetricLoading(false);
  };

  const addMetricView = (
    <select
      defaultValue="add-metric"
      defaultChecked={true}
      onChange={(e) => {
        if (e.target.value === "atom-atom") {
          newAtomAtomDistance();
        }
        // Re-select the default value
        e.target.value = "add-metric";
      }}
      value={"add-metric"}
    >
      <option key="add-metric" value={"add-metric"} disabled>
        Add metric
      </option>
      <option key="new-metric-atom" value={"atom-atom"}>
        Atom-Atom distance
      </option>
    </select>
  );

  return (
    <div className="blurred-squared-div">
      <div
        style={{
          display: "flex",
          justifyContent: "space-between",
          alignItems: "center",
        }}
      >
        <div className="title-color">Compare structures</div>
        <div
          style={{
            display: "flex",
            gap: "1rem",
            alignItems: "center",
          }}
        >
          {metricLoading ? <div className="loader"></div> : addMetricView}
          <div>
            <button
              onClick={() => {
                // Reset columns
                gridRef.current?.columnApi.resetColumnState();
              }}
            >
              Reset columns
            </button>
          </div>
        </div>
      </div>
      {/* On div wrapping Grid a) specify theme CSS Class Class and b) sets Grid size */}
      <div className="ag-theme-alpine">
        <AgGridReact
          ref={gridRef} // Ref for accessing Grid's API
          className="pele-table" // CSS Styles for AgGridReact
          rowData={rowData} // Row Data for Rows
          columnDefs={columnDefs} // Column Defs for Columns
          defaultColDef={defaultColDef} // Default Column Properties
          animateRows={true} // Optional - set to 'true' to have rows animate when sorted
          rowSelection="multiple" // Options - allows click selection of rows
          onCellClicked={cellClickedListener} // Optional - registering for Grid Event
          getRowId={(params) => params.data.id}
          domLayout="autoHeight"
        />
      </div>
    </div>
  );
}

interface PELEInputProps {
  nbdData: NBDSuiteData;
}

function PELEInput(props: PELEInputProps) {
  const [info, setInfo] = useState<any>(null);

  async function loadInputInfo() {
    await props.nbdData.getInputInfo();
    setInfo(props.nbdData.inputInfo);
  }

  useEffect(() => {
    if (props.nbdData) {
      loadInputInfo();
    }
  }, []);

  return (
    <div>
      <div className="blurred-squared-div">
        <div className="title-color">Input file</div>
        <textarea
          value={info || ""}
          onChange={(e) => setInfo(e.target.value)}
          style={{ minHeight: 500, height: "100%", width: "100%" }}
        />
      </div>
    </div>
  );
}

function NBDSuiteResults() {
  const [openedFile, setOpenedFile] = useState("");

  const [complexes, setComplexes] = useState<Array<string>>([]);
  const [data, setData] = useState([]);

  const [selectedComplex, setSelectedComplex] = useState("");

  const [error, setError] = useState(false);

  const [extraData, setExtraData] = useState<Array<any>>([]);

  const errorMsg = useRef("");
  const NBDData = useRef<NBDSuiteData | null>(null);

  const addExtraData = (data: any) => {
    setExtraData({ ...extraData, ...data });
  };

  const updateSelectedComplex = (e: any) => {
    if (NBDData.current) {
      NBDData.current!.selectedComplex = e.target.value;
    }
    setSelectedComplex(e.target.value);
  };

  const fetchNBDData = async () => {
    const nbddata = new NBDSuiteData(openedFile);

    try {
      await nbddata.getSimulationName();

      await nbddata.getComplexes();

      setComplexes(nbddata.complexes);

      await nbddata.getPlotData(nbddata.complexes[0]);

      setData(nbddata.plotData);

      nbddata.getInputPDB();

      NBDData.current = nbddata;
    } catch (e) {
      setError(true);
      errorMsg.current = e.message;
      return;
    }
  };

  useEffect(() => {
    if (openedFile) {
      fetchNBDData();
    }
  }, [openedFile]);

  return (
    <div>
      {openedFile && (
        <div>
          {error && <p>Error reading results {errorMsg.current}</p>}
          {NBDData.current && (
            <div className="flex-center-col">
              <PELEPlot
                nbdData={NBDData.current}
                addExtraData={addExtraData}
                updateSelectedComplex={updateSelectedComplex}
              />
              <PELETable nbdData={NBDData.current} extraData={extraData} />
              <PELEInput nbdData={NBDData.current!} />
            </div>
          )}
        </div>
      )}
      <div
        style={{
          height: openedFile ? "100%" : "100vh",
          marginBottom: openedFile ? "2rem" : 0,
        }}
      >
        <OpenInputFile setOpenedFile={setOpenedFile} />
      </div>
    </div>
  );
}

export { NBDSuiteResults };
