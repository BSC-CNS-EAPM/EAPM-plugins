const peleSimulationFolder = document.querySelector('#peleSimulationFolder')
const peleSimulationFolderBtn = document.querySelector('#peleSimulationFolderBtn')
const peleOutputFolder = document.querySelector('#peleOutputFolder')
const peleForm = document.querySelector('#peleForm')
const plotForm = document.querySelector('#plotForm')
const proteinSelector = document.querySelector('#proteinSelector')
const ligandSelector = document.querySelector('#ligandSelector')
const distanceSelector = document.querySelector('#distanceSelector')
const metricSelector = document.querySelector('#metricSelector')
const colorSelector = document.querySelector('#colorSelector')
const loader = document.querySelector('#loader')

const byMetric = document.querySelector('#byMetric')
const colorByMetric = document.querySelector('#colorByMetric')
const filterByMetric = document.querySelector('#filterByMetric')
const metricSerL = document.querySelector('#metric_SER-L-slider')
const metricSerHis = document.querySelector('#metric_SER-HIS-slider')
const metricHisAsp = document.querySelector('#metric_HIS-ASP-slider')
const metricSerLOutput = document.querySelector('#metric_SER-L-output')
const metricSerHisOutput = document.querySelector('#metric_SER-HIS-output')
const metricHisAspOutput = document.querySelector('#metric_HIS-ASP-output')
const metricSerLLabel = document.querySelector('#metric_SER-L-label')
const metricSerHisLabel = document.querySelector('#metric_SER-HIS-label')
const metricHisAspLabel = document.querySelector('#metric_HIS-ASP-label')

let fetchData
let dataColumn
let gridApi

metricSerL.addEventListener('change', () => {
    metricSerLOutput.value = metricSerL.value
})

metricSerHis.addEventListener('change', () => {
    metricSerHisOutput.value = metricSerHis.value
})

metricHisAsp.addEventListener('change', () => {
    metricHisAspOutput.value = metricHisAsp.value
})


const horusData = window.parent.extensionData;

if(horusData?.peleFolder !== undefined){
    peleSimulationFolder.innerHTML = horusData.peleFolder
}
// TODO pendiente refactorizar para no tener que crear el evento de plotly_click


// parent.molstar.loadPDB(pdbData, label)

const addOption = (selector, array) => {
    array.forEach(i => {
        option = document.createElement('option')
        option.id = i
        option.value = i
        option.innerHTML = i
        selector.append(option)
    });
}

const addDictOption = (selector, array) => {
    Object.keys(array).forEach(i => {
        option = document.createElement('option')
        option.id = i
        option.value = i
        option.innerHTML = i
        selector.append(option)
    });
}

peleSimulationFolderBtn.addEventListener("click", (e) => {
    e.preventDefault()

    const href = window.location.href

    const data = {
        peleSimulationFolder: peleSimulationFolder.value,
        peleOutputFolder: peleOutputFolder.value,
        desiredProtein: null,
        desiredLigand: null
    }

    loader.style.display = 'block'
    fetch(href + '/peleResults', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify(data)
    })
        .then(response => response.json())
        .then(parsedData => {
            if(parsedData.ok){
                fetchData = parsedData
                console.log(parsedData)

                const proteins = parsedData.distances
                const proteinValue = Object.keys(proteins)[0]
                const ligands = parsedData.distances[proteinValue]
                const distances = parsedData.distances[proteinValue][Object.keys(ligands)[0]]
                const metrics = parsedData.metrics.map(metric => "metric_" + metric);

                addDictOption(proteinSelector, proteins)
                addDictOption(ligandSelector, ligands)
                addOption(distanceSelector, distances)
                addOption(metricSelector, metrics)

                peleForm.style.display = 'none'
                plotForm.style.display = 'block'
                loader.style.display = 'none'

                showPlot(distanceSelector.value, "Binding Energy")
                createGrid()

                pelePlot.on("plotly_click", (data) => {
                    console.log(gridApi.rowModel.rowsToDisplay)

                    let newData = []

                    gridApi.rowModel.rowsToDisplay.forEach(row => {
                        newData.push(row.data)
                    });

                    const pointData = data.points[0].data.data;

                    newData.push({
                        "Protein": proteinSelector.value,
                        "Ligand": ligandSelector.value,
                        "Epoch": pointData["Epoch"],
                        "Trajectory": pointData["Trajectory"],
                        "Accepted Pele Steps": pointData["Accepted Pele Steps"],
                        "Step": pointData["Step"],
                        "Total Energy": pointData["Total Energy"],
                        "Binding Energy": pointData["Binding Energy"],
                        "Ligand SASA": pointData["Ligand SASA"],
                        "metric_SER-L": pointData["metric_SER-L"],
                        "metric_SER-HIS": pointData["metric_SER-HIS"],
                        "metric_HIS-ASP": pointData["metric_HIS-ASP"]
                    })

                    gridApi.setGridOption('rowData', newData)
                });
            }else{
                console.log('Error: ' + parsedData.msg)
            }
        })
        .catch(error => {
            console.error('Error en la solicitud:', error);
            console.log('JSON del error:', error.response);
        });
    
})

const createGrid = () => {
    const columnNames = [
        "Protein",
        "Ligand",
        "Epoch",
        "Trajectory",
        "Accepted Pele Steps",
        "Step",
        "Total Energy",
        "Binding Energy",
        "Ligand SASA",
        "metric_SER-L",
        "metric_SER-HIS",
        "metric_HIS-ASP"
    ];

    const columnDefs = columnNames.map(column => {
        return { field: column };
    });
    
    const gridOptions = {
        rowData: [],
        columnDefs: columnDefs,
        rowSelection: 'single',
        onRowClicked: event => {
            console.log(event.data);
            const href = window.location.href;
        
            const jsonData = JSON.stringify(event.data);
        
            fetch(href + '/plotlyClick', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: jsonData
            })
            .then(response => response.json())
            .then(parsedData => {
                if (parsedData.ok) {
                    console.log(parsedData.msg);
                }
            });
        }
    };
    
    const gridElement = document.querySelector('#myGrid');
    gridApi = agGrid.createGrid(gridElement, gridOptions);
}

const changeDisplayStyle = (array, string) => {
    array.forEach(element => {
        element.style.display = string
    });
}

const showPlot = (xVal, yVal) => {
    let toPlotData = []

    for (const value in fetchData.dict) {
        bindingEnergy = fetchData.dict[value]["Binding Energy"]
        acceptedPeleStep = fetchData.dict[value]["Accepted Pele Steps"]

        let marker

        if(filterByMetric.checked){
            if((fetchData.dict[value]["metric_SER-L"] > metricSerL.value) 
            || (fetchData.dict[value]["metric_SER-HIS"] > metricSerHis.value) 
            || (fetchData.dict[value]["metric_HIS-ASP"] > metricHisAsp.value)){
                continue
            }
        }

        if(colorSelector.value == 'None'){
            marker = {
                size: 6,
                color: colorByMetric.checked ? (fetchData.dict[value]["metric_SER-L"] <= metricSerL.value) 
                && (fetchData.dict[value]["metric_SER-HIS"] <= metricSerHis.value) 
                && (fetchData.dict[value]["metric_HIS-ASP"] <= metricHisAsp.value) 
                ? 'rgb(255, 0, 0)' : 'rgb(0, 0, 0)' : 'rgba(45, 85, 255, 0.7)'
            }
        }else{
            marker = {
                size: 6,
                color: colorSelector.value,
                colorscale: 'Viridis',
                colorbar: {
                    title: colorSelector.value
                }
             }
        }

        const data = {
            x: [fetchData.dict[value][xVal]],
            y: [fetchData.dict[value][yVal]],
            type: "scatter",
            mode: "markers",
            marker: marker,
            data: fetchData.dict[value],
            showlegend: false
        }

        toPlotData.push(data)
    }

    console.log(toPlotData)

    var pelePlot = document.getElementById('pelePlot'),
    layout = {
        title: 'Scatter plot',
        hovermode: 'closest',
        width: 0.60 * pelePlot.offsetWidth,
        height: 0.60 * pelePlot.offsetWidth,
        xaxis: {
            title: xVal
        },
        yaxis: {
            title: yVal
        }
    };

    Plotly.newPlot('pelePlot', toPlotData, layout);
}

plotForm.addEventListener("change", (e) => {

    if(e.target.id === "proteinSelector"){

        const href = window.location.href

        const data = {
            peleSimulationFolder: peleSimulationFolder.value,
            peleOutputFolder: peleOutputFolder.value,
            desiredProtein: proteinSelector.value,
            desiredLigand: Object.keys(fetchData.distances[proteinSelector.value])[0]
        }

        proteinSelector.disabled = true
        ligandSelector.disabled = true
        distanceSelector.disabled = true
        colorSelector.disabled = true
        loader.style.display = 'block'
        fetch(href + '/peleResults', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(data)
        })
            .then(response => response.json())
            .then(parsedData => {
                if(parsedData.ok){
                    fetchData = parsedData
                    console.log(parsedData)
                    
                    ligandSelector.innerHTML = ''
                    distanceSelector.innerHTML = ''

                    const proteinValue = e.target.value
                    const ligands = parsedData.distances[proteinValue]
                    
                    addDictOption(ligandSelector, ligands)
                    
                    const distances = parsedData.distances[proteinValue][ligandSelector.value]

                    addOption(distanceSelector, distances)

                    proteinSelector.disabled = false
                    ligandSelector.disabled = false
                    distanceSelector.disabled = false
                    colorSelector.disabled = false
                    loader.style.display = 'none'

                    showPlot(distanceSelector.value, "Binding Energy")

                    pelePlot.on("plotly_click", (data) => {
                        console.log(gridApi.rowModel.rowsToDisplay)
    
                        let newData = []
    
                        gridApi.rowModel.rowsToDisplay.forEach(row => {
                            newData.push(row.data)
                        });
    
                        const pointData = data.points[0].data.data;
    
                        newData.push({
                            "Protein": proteinSelector.value,
                            "Ligand": ligandSelector.value,
                            "Epoch": pointData["Epoch"],
                            "Trajectory": pointData["Trajectory"],
                            "Accepted Pele Steps": pointData["Accepted Pele Steps"],
                            "Step": pointData["Step"],
                            "Total Energy": pointData["Total Energy"],
                            "Binding Energy": pointData["Binding Energy"],
                            "Ligand SASA": pointData["Ligand SASA"],
                            "metric_SER-L": pointData["metric_SER-L"],
                            "metric_SER-HIS": pointData["metric_SER-HIS"],
                            "metric_HIS-ASP": pointData["metric_HIS-ASP"]
                        })
    
                        gridApi.setGridOption('rowData', newData)
                    });
                }else{
                    console.log('Error: ' + parsedData.msg)
                }
            })
            .catch(error => console.error('Error en la solicitud:', error));
    }else if(e.target.id === "ligandSelector"){
        const href = window.location.href

        const data = {
            peleSimulationFolder: peleSimulationFolder.value,
            peleOutputFolder: peleOutputFolder.value,
            desiredProtein: proteinSelector.value,
            desiredLigand: ligandSelector.value
        }

        proteinSelector.disabled = true
        ligandSelector.disabled = true
        distanceSelector.disabled = true
        colorSelector.disabled = true
        loader.style.display = 'block'
        fetch(href + '/testEndpoint', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(data)
        })
            .then(response => response.json())
            .then(parsedData => {
                if(parsedData.ok){
                    fetchData = parsedData
                    console.log(parsedData)

                    const proteinValue = proteinSelector.value
                    const distances = fetchData.distances[proteinValue][ligandSelector.value]

                    distanceSelector.innerHTML = ''

                    addOption(distanceSelector, distances)

                    proteinSelector.disabled = false
                    ligandSelector.disabled = false
                    distanceSelector.disabled = false
                    colorSelector.disabled = false
                    loader.style.display = 'none'

                    showPlot(distanceSelector.value, "Binding Energy")

                    pelePlot.on("plotly_click", (data) => {
                        console.log(gridApi.rowModel.rowsToDisplay)
    
                        let newData = []
    
                        gridApi.rowModel.rowsToDisplay.forEach(row => {
                            newData.push(row.data)
                        });
    
                        const pointData = data.points[0].data.data;
    
                        newData.push({
                            "Protein": proteinSelector.value,
                            "Ligand": ligandSelector.value,
                            "Epoch": pointData["Epoch"],
                            "Trajectory": pointData["Trajectory"],
                            "Accepted Pele Steps": pointData["Accepted Pele Steps"],
                            "Step": pointData["Step"],
                            "Total Energy": pointData["Total Energy"],
                            "Binding Energy": pointData["Binding Energy"],
                            "Ligand SASA": pointData["Ligand SASA"],
                            "metric_SER-L": pointData["metric_SER-L"],
                            "metric_SER-HIS": pointData["metric_SER-HIS"],
                            "metric_HIS-ASP": pointData["metric_HIS-ASP"]
                        })
    
                        gridApi.setGridOption('rowData', newData)
                    });
                }else{
                    console.log('Error: ' + parsedData.msg)
                }
            })
            .catch(error => console.error('Error en la solicitud:', error));
    }else if(e.target.id === "distanceSelector"){
        showPlot(distanceSelector.value, "Binding Energy")
        pelePlot.on("plotly_click", (data) => {
            console.log(gridApi.rowModel.rowsToDisplay)

            let newData = []

            gridApi.rowModel.rowsToDisplay.forEach(row => {
                newData.push(row.data)
            });

            const pointData = data.points[0].data.data;

            newData.push({
                "Protein": proteinSelector.value,
                "Ligand": ligandSelector.value,
                "Epoch": pointData["Epoch"],
                "Trajectory": pointData["Trajectory"],
                "Accepted Pele Steps": pointData["Accepted Pele Steps"],
                "Step": pointData["Step"],
                "Total Energy": pointData["Total Energy"],
                "Binding Energy": pointData["Binding Energy"],
                "Ligand SASA": pointData["Ligand SASA"],
                "metric_SER-L": pointData["metric_SER-L"],
                "metric_SER-HIS": pointData["metric_SER-HIS"],
                "metric_HIS-ASP": pointData["metric_HIS-ASP"]
            })

            gridApi.setGridOption('rowData', newData)
        });
    }else if(e.target.id === 'colorByMetric'){

        let metricsDisplay = [
            metricSerL, metricSerLOutput, metricSerLLabel, metricSerHis, metricSerHisOutput, 
            metricSerHisLabel, metricHisAsp, metricHisAspOutput, metricHisAspLabel
        ] 

        if(e.target.checked){
            changeDisplayStyle(metricsDisplay, 'inline-block')
        }else{
            changeDisplayStyle(metricsDisplay, 'none')
        }

        showPlot(distanceSelector.value, "Binding Energy")
        pelePlot.on("plotly_click", (data) => {
            console.log(gridApi.rowModel.rowsToDisplay)

            let newData = []

            gridApi.rowModel.rowsToDisplay.forEach(row => {
                newData.push(row.data)
            });

            const pointData = data.points[0].data.data;

            newData.push({
                "Protein": proteinSelector.value,
                "Ligand": ligandSelector.value,
                "Epoch": pointData["Epoch"],
                "Trajectory": pointData["Trajectory"],
                "Accepted Pele Steps": pointData["Accepted Pele Steps"],
                "Step": pointData["Step"],
                "Total Energy": pointData["Total Energy"],
                "Binding Energy": pointData["Binding Energy"],
                "Ligand SASA": pointData["Ligand SASA"],
                "metric_SER-L": pointData["metric_SER-L"],
                "metric_SER-HIS": pointData["metric_SER-HIS"],
                "metric_HIS-ASP": pointData["metric_HIS-ASP"]
            })

            gridApi.setGridOption('rowData', newData)
        })
    }else if(e.target.id === 'metric_SER-L-slider' || e.target.id === 'metric_SER-HIS-slider' || e.target.id === 'metric_HIS-ASP-slider'){
        showPlot(distanceSelector.value, "Binding Energy")
        pelePlot.on("plotly_click", (data) => {
            console.log(gridApi.rowModel.rowsToDisplay)

            let newData = []

            gridApi.rowModel.rowsToDisplay.forEach(row => {
                newData.push(row.data)
            });

            const pointData = data.points[0].data.data;

            newData.push({
                "Protein": proteinSelector.value,
                "Ligand": ligandSelector.value,
                "Epoch": pointData["Epoch"],
                "Trajectory": pointData["Trajectory"],
                "Accepted Pele Steps": pointData["Accepted Pele Steps"],
                "Step": pointData["Step"],
                "Total Energy": pointData["Total Energy"],
                "Binding Energy": pointData["Binding Energy"],
                "Ligand SASA": pointData["Ligand SASA"],
                "metric_SER-L": pointData["metric_SER-L"],
                "metric_SER-HIS": pointData["metric_SER-HIS"],
                "metric_HIS-ASP": pointData["metric_HIS-ASP"]
            })

            gridApi.setGridOption('rowData', newData)
        })
    }else if(e.target.id === 'colorSelector'){
        showPlot(distanceSelector.value, "Binding Energy")
        pelePlot.on("plotly_click", (data) => {
            console.log(gridApi.rowModel.rowsToDisplay)

            let newData = []

            gridApi.rowModel.rowsToDisplay.forEach(row => {
                newData.push(row.data)
            });

            const pointData = data.points[0].data.data;

            newData.push({
                "Protein": proteinSelector.value,
                "Ligand": ligandSelector.value,
                "Epoch": pointData["Epoch"],
                "Trajectory": pointData["Trajectory"],
                "Accepted Pele Steps": pointData["Accepted Pele Steps"],
                "Step": pointData["Step"],
                "Total Energy": pointData["Total Energy"],
                "Binding Energy": pointData["Binding Energy"],
                "Ligand SASA": pointData["Ligand SASA"],
                "metric_SER-L": pointData["metric_SER-L"],
                "metric_SER-HIS": pointData["metric_SER-HIS"],
                "metric_HIS-ASP": pointData["metric_HIS-ASP"]
            })

            gridApi.setGridOption('rowData', newData)
        })
    }else if(e.target.id === 'filterByMetric'){

        let metricsDisplay = [
            metricSerL, metricSerLOutput, metricSerLLabel, metricSerHis, metricSerHisOutput, 
            metricSerHisLabel, metricHisAsp, metricHisAspOutput, metricHisAspLabel
        ] 

        if(e.target.checked){
            changeDisplayStyle(metricsDisplay, 'inline-block')
        }else{
            changeDisplayStyle(metricsDisplay, 'none')
        }

        showPlot(distanceSelector.value, "Binding Energy")
        pelePlot.on("plotly_click", (data) => {
            console.log(gridApi.rowModel.rowsToDisplay)

            let newData = []

            gridApi.rowModel.rowsToDisplay.forEach(row => {
                newData.push(row.data)
            });

            const pointData = data.points[0].data.data;

            newData.push({
                "Protein": proteinSelector.value,
                "Ligand": ligandSelector.value,
                "Epoch": pointData["Epoch"],
                "Trajectory": pointData["Trajectory"],
                "Accepted Pele Steps": pointData["Accepted Pele Steps"],
                "Step": pointData["Step"],
                "Total Energy": pointData["Total Energy"],
                "Binding Energy": pointData["Binding Energy"],
                "Ligand SASA": pointData["Ligand SASA"],
                "metric_SER-L": pointData["metric_SER-L"],
                "metric_SER-HIS": pointData["metric_SER-HIS"],
                "metric_HIS-ASP": pointData["metric_HIS-ASP"]
            })

            gridApi.setGridOption('rowData', newData)
        })
    }else if(e.target.id === 'byMetric'){
        if(e.target.checked){
            changeDisplayStyle([distanceSelector], 'none')
            changeDisplayStyle([metricSelector], 'inline-block')
        }else{
            changeDisplayStyle([distanceSelector], 'inline-block')
            changeDisplayStyle([metricSelector], 'none')
        }
    }else if(e.target.id === 'metricSelector'){
        showPlot(metricSelector.value, "Binding Energy")
        pelePlot.on("plotly_click", (data) => {
            console.log(gridApi.rowModel.rowsToDisplay)

            let newData = []

            gridApi.rowModel.rowsToDisplay.forEach(row => {
                newData.push(row.data)
            });

            const pointData = data.points[0].data.data;

            newData.push({
                "Protein": proteinSelector.value,
                "Ligand": ligandSelector.value,
                "Epoch": pointData["Epoch"],
                "Trajectory": pointData["Trajectory"],
                "Accepted Pele Steps": pointData["Accepted Pele Steps"],
                "Step": pointData["Step"],
                "Total Energy": pointData["Total Energy"],
                "Binding Energy": pointData["Binding Energy"],
                "Ligand SASA": pointData["Ligand SASA"],
                "metric_SER-L": pointData["metric_SER-L"],
                "metric_SER-HIS": pointData["metric_SER-HIS"],
                "metric_HIS-ASP": pointData["metric_HIS-ASP"]
            })

            gridApi.setGridOption('rowData', newData)
        })
    }
})


