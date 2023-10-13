import { horusPost } from "../../../../Views/Utils/utils";
import { AddButton, DeleteButton, SaveButton } from "./renderers";
import downloadjs from "downloadjs";

class NBDSuiteData {
  /**
   * @description - path to the NBDSuite folder
   * @type {string}
   */
  path: string;

  /**
   * @description - The list of data to plot
   * @type {any}
   */
  plotData: any;

  /**
   * @description - The list of complexes in the opened NBDSuite folder
   * @type {Array<string>}
   */
  complexes: Array<string>;

  /**
   * @description - The name of the opened NBDSuite simulation
   * @type {any}
   * @memberof NBDSuiteData
   */
  inputInfo: any;

  /**
   * @description - The selected complex to plot
   * @type {string}
   * @memberof NBDSuiteData
   * */
  selectedComplex: string;

  /**
   * @description - The name of the simulation
   * @type {string}
   * @memberof NBDSuiteData
   * */
  simulationName: string;

  /**
   *
   * @param path - path to the NBDSuite folder
   */
  constructor(path: string) {
    this.path = path;
  }

  async getSimulationName() {
    const body = JSON.stringify({ path: this.path });

    const header = {
      "Content-Type": "application/json",
      Accept: "application/json",
    };

    const href = window.location.href;

    const postTo = href + "getInputSimulationName";

    const request = await horusPost(postTo, header, body);

    const data = await request.json();

    if (!data.ok) {
      throw new Error("Error getting input info");
    }

    this.simulationName = data.data;
  }

  async getInputInfo() {
    const body = JSON.stringify({ path: this.path });

    const header = {
      "Content-Type": "application/json",
      Accept: "application/json",
    };

    const href = window.location.href;

    const postTo = href + "getInputInfo";

    const request = await horusPost(postTo, header, body);

    const data = await request.json();

    if (!data.ok) {
      throw new Error("Error getting input info");
    }

    this.inputInfo = data.data;
  }

  async getComplexes() {
    const body = JSON.stringify({ path: this.path });

    const header = {
      "Content-Type": "application/json",
      Accept: "application/json",
    };

    const postTo = window.location.href + "loadComplexes";

    const request = await horusPost(postTo, header, body);

    const data = await request.json();

    // Check the response data
    if (!data.ok) {
      throw new Error("Error getting plot data");
    }

    this.complexes = data.complexes;

    this.selectedComplex = this.complexes[0];
  }

  async getPlotData(complex: string) {
    // The POST request url
    const url = window.location.href + "getPlotData";

    // The POST request body
    const body = JSON.stringify({
      path: this.path,
      complex: complex,
    });

    // The POST request header
    const header = {
      "Content-Type": "application/json",
      Accept: "application/json",
    };

    // The POST request options
    const response = await horusPost(url, header, body);

    // Check the response status
    if (response.status != 200) {
      throw new Error("Error getting plot data");
    }

    // Get the response data
    const data = await response.json();

    // Check the response data
    if (!data.ok) {
      throw new Error(data.msg);
    }

    // Set the plot data
    this.plotData = data.plotdata;
  }

  /**
   * @description - Returns the possible x/y-axis options based on the plot data
   * @returns {Array<string>}
   * @memberof NBDSuiteData
   */
  axisOptions() {
    const options: Array<string> = [];

    // Get one of the plot data objects
    const someData = this.plotData[Object.keys(this.plotData)[0]][0];

    if (!someData) {
      throw new Error("No plot data");
    }

    // Loop through the plot data
    for (const key in someData) {
      if (
        key == "id" ||
        key == "filename" ||
        key == "filepath" ||
        key == "trajectory"
      ) {
        continue;
      }
      // Add the key to the options
      options.push(key);
    }

    return options;
  }

  /**
   * @description - Returns the string of the PDB file used in the input
   * @returns {string}
   * @memberof NBDSuiteData
   */
  async getInputPDB() {
    const body = JSON.stringify({ path: this.path });

    const header = {
      "Content-Type": "application/json",
      Accept: "application/json",
    };

    const href = window.location.href;

    const postTo = href + "loadInputPDB";

    const request = await horusPost(postTo, header, body);

    const data = await request.json();

    if (!data.ok) {
      throw new Error("Error getting input PDB");
    }

    const pdb = data.data.pdb;
    const name = data.data.name;

    // Retrive the molstar object from the window
    window.parent.molstar?.loadPDBString(pdb, name);
  }

  /**
   * @description - Returns the plot of the data
   * @returns {string}
   * @memberof NBDSuiteData
   **/
  getPlot(xAxis: string, yAxis: string) {
    const clusterColors = {
      A: "#636EFA",
      B: "#EF553B",
      C: "#00CC96",
      D: "#AB63FA",
      E: "#FFA15A",
      Other: "#909091",
    };

    const dToP: Array<any> = [];
    const legendEntries: any = {};
    for (const d in this.plotData) {
      const currentData = this.plotData[d];
      for (const cd in currentData) {
        const toplot = currentData[cd];

        const xVal = toplot[xAxis];
        const yVal = toplot[yAxis];
        const cluster = toplot["Cluster label"];
        let symbol = "circle";

        if (d === "repr") {
          symbol = "star";
        }

        let showLegend = false;
        if (!legendEntries[cluster]) {
          legendEntries[cluster] = {
            color: clusterColors[cluster] || "red",
            symbol: symbol,
          };
          showLegend = true;
        }

        const data = {
          x: [xVal],
          y: [yVal],
          type: "scatter",
          mode: "markers",
          marker: {
            color: clusterColors[cluster] || clusterColors.Other,
            symbol: symbol,
          },
          name: cluster,
          // Hide the legend for the repr data
          showlegend: showLegend,
          data: toplot,
        };

        dToP.push(data);
      }
    }
    return dToP;
  }

  /**
   * @description - Returns the initial table data. (Representative data)
   * @returns {Array<any>}
   * @memberof NBDSuiteData
   **/
  getTable() {
    const columns = this.axisOptions().map((c: string) => {
      if (
        c === "trajectory" ||
        c === "id" ||
        c === "filename" ||
        c === "filepath" ||
        c === "Type"
      ) {
        // Do not show the trajectory column
        return { headerName: c.toUpperCase(), field: c, hide: true };
      }
      return { headerName: c, field: c, resizable: true, flex: 1 };
    });

    // Add a "View" column
    columns.unshift({
      headerName: "Add to view",
      field: "show",
      cellRenderer: AddButton,
      resizable: true,
      flex: 1,
    });

    // Add a download/save column
    columns.push({
      headerName: "Save",
      field: "save",
      flex: 1,
      cellRenderer: SaveButton,
    });

    // Add a delete column
    columns.push({
      headerName: "Delete",
      field: "delete",
      flex: 1,
      cellRenderer: DeleteButton,
    });

    // Move the "Cluster label" column to the 2nd position
    const clusterLabelIndex = columns.findIndex(
      (c: any) => c.field === "Cluster label"
    );
    const clusterLabel = columns.splice(clusterLabelIndex, 1)[0];

    columns.splice(1, 0, clusterLabel);

    const rows = this.plotData.repr;
    // Add to the repr data a delete column
    rows.forEach((d: any) => {
      d.delete = "Delete";
    });

    return { columns: columns, rows: rows };
  }

  private async retrivePDBContentsFromRow(row: any) {
    const filepath = row["filepath"];
    const filename = row["filename"];
    const trajectory = row["trajectory"];
    const step = row["Step"];

    if (!filepath) {
      alert("PDBs of non-representative data cannot be viewed yet");
      return;
    }

    const header = {
      "Content-Type": "application/json",
      Accept: "application/json",
    };

    const body = JSON.stringify({
      path: this.path,
      selectedComplex: this.selectedComplex,
      filepath: filepath,
      filename: filename,
      trajectory: trajectory,
      step: step,
    });

    const href = window.location.href;

    const postTo = href + "getPDB";

    // Get the pdb file
    const response = await horusPost(postTo, header, body);

    // Get the response data
    const data = await response.json();

    // Check the response data
    if (!data.ok) {
      alert(data.msg);
    }

    // Get the pdb data (pdb and name)
    return data.data;
  }

  /**
   * @description - Fetches the given PDB file and adds it to molstar.
   *  Untruncates the file if necessary
   * @param {string} filepath
   * @param {string} name
   * @memberof NBDSuiteData
   * @returns {void}
   */
  async addPDBToMolstarFromRow(row: any) {
    const filename = row["filename"];
    const pdb = await this.retrivePDBContentsFromRow(row);

    if (!pdb) {
      return;
    }

    // Add the pdb to molstar
    window.parent.molstar?.loadPDBString(pdb.pdb, filename);
  }

  async savePDBfromRow(row: any) {
    const pdb = await this.retrivePDBContentsFromRow(row);

    if (!pdb) {
      return;
    }

    const header = {
      "Content-Type": "application/json",
      Accept: "application/json",
    };

    const body = JSON.stringify({
      contents: pdb.pdb,
      filename: pdb.name,
    });

    const postTo = "/savecontents";

    // Check if the response can be parsed as JSON (ok file saved on desktop) or as a blob (we are on server mode)
    try {
      // Get the pdb file
      const response = await horusPost(postTo, header, body);
      const data = await response.json();
      if (!data.ok) {
        alert(data.msg);
      }
      return;
    } catch (e) {
      try {
        // If the POST request fails, try to download the file as a blob
        const newHeader = {
          "Content-Type": "application/octet-stream",
          Accept: "application/octet-stream",
        };

        const response = await horusPost(postTo, header, body);

        const blob = await response.blob();

        downloadjs(blob, pdb.name, "text/plain");
      } catch (e) {
        alert("Unable to save the file");
      }
    }
  }

  async addAtomAtomDistance() {
    // Get from mosltar the selected atoms

    const selectedAtoms = window.parent.molstar?.getSelectedStructures();

    // Send the selected atoms to the server
    const header = {
      "Content-Type": "application/json",
      Accept: "application/json",
    };

    const body = JSON.stringify({
      selectedAtoms: selectedAtoms,
      selectedComplex: this.selectedComplex,
    });

    const postTo = window.location.href + "getAtomAtomDistance";

    const response = await horusPost(postTo, header, body);

    const data = await response.json();

    if (!data.ok) {
      alert(data.msg);
    }

    const atomAtomDistance = data.data;

    // Add the atom-atom distance to the plot data
    // ** WIP **
  }
}

export default NBDSuiteData;
