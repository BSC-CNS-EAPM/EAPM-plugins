const peleSimulationFolder = document.querySelector('#peleSimulationFolder')
const peleSimulationFolderBtn = document.querySelector('#peleSimulationFolderBtn')
const peleOutputFolder = document.querySelector('#peleOutputFolder')
const peleForm = document.querySelector('#peleForm')
const plotForm = document.querySelector('#plotForm')
const proteinSelector = document.querySelector('#proteinSelector')
const ligandSelector = document.querySelector('#ligandSelector')
const distanceSelector = document.querySelector('#distanceSelector')
const colorSelector = document.querySelector('#colorSelector')
const loader = document.querySelector('#loader')




let fetchData
let dataColumn

const horusData = window.parent.extensionData;

if(horusData?.peleFolder !== undefined){
    peleSimulationFolder.innerHTML = horusData.peleFolder
}

//   parent.molstar.loadPDB(pdbData, label)

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
                
                addDictOption(proteinSelector, proteins)

                addDictOption(ligandSelector, ligands)

                addOption(distanceSelector, distances)

                console.log(distances)

                peleForm.style.display = 'none'

                plotForm.style.display = 'block'

                const bindingEnergyIndex = fetchData.dict.columns.indexOf("Binding Energy")

                dataColumn = fetchData.dict.data.map(row => row[bindingEnergyIndex])

                loader.style.display = 'none'

                showPlot()

                pelePlot.on("plotly_click", (data) => {
                    const d = data.points[0];

                    const distance = d.x 
                    const bindingEnergy = d.y

                    const href = window.location.href

                    const plotlyData = {
                        distance: distance,
                        bindingEnergy: bindingEnergy
                    }

                    fetch(href + '/plotlyClick', {
                        method: 'POST',
                        headers: {
                            'Content-Type': 'application/json'
                        },
                        body: JSON.stringify(plotlyData)
                    })
                        .then(response => response.json())
                        .then(parsedData => {
                            if(parsedData.ok){
                                console.log(parsedData.msg)
                            }
                        })

                    console.log(d)
                    // addRowToTable(d);
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

const showGrid = () => {
    const proteinColumnIndex = fetchData.dict.columns.indexOf('Protein');
    const proteinColumn = fetchData.dict.data.map(row => row[proteinColumnIndex]);

    const ligandColumnIndex = fetchData.dict.columns.indexOf('Ligand');
    const ligandColumn = fetchData.dict.data.map(row => row[ligandColumnIndex]);
    
    const epochColumnIndex = fetchData.dict.columns.indexOf('Epoch');
    const epochColumn = fetchData.dict.data.map(row => row[epochColumnIndex]);

    const trajectoryColumnIndex = fetchData.dict.columns.indexOf('Trajectory');
    const trajectoryColumn = fetchData.dict.data.map(row => row[trajectoryColumnIndex]);

    const acceptedPeleStepsColumnIndex = fetchData.dict.columns.indexOf('Accepted Pele Steps');
    const acceptedPeleStepsColumn = fetchData.dict.data.map(row => row[acceptedPeleStepsColumnIndex]);

    const stepColumnIndex = fetchData.dict.columns.indexOf('Step');
    const stepColumn = fetchData.dict.data.map(row => row[stepColumnIndex]);

    const totalEnergyColumnIndex = fetchData.dict.columns.indexOf('Total Energy');
    const totalEnergyColumn = fetchData.dict.data.map(row => row[totalEnergyColumnIndex]);

    const bindingEnergyColumnIndex = fetchData.dict.columns.indexOf('Binding Energy');
    const bindingEnergyColumn = fetchData.dict.data.map(row => row[bindingEnergyColumnIndex]);

    const ligandSASAColumnIndex = fetchData.dict.columns.indexOf('Ligand SASA');
    const ligandSASAColumn = fetchData.dict.data.map(row => row[ligandSASAColumnIndex]);

    const metricSERLColumnIndex = fetchData.dict.columns.indexOf('metric_SER-L');
    const metricSERLColumn = fetchData.dict.data.map(row => row[metricSERLColumnIndex]);

    const metricSERHISColumnIndex = fetchData.dict.columns.indexOf('metric_SER-HIS');
    const metricSERHISColumn = fetchData.dict.data.map(row => row[metricSERHISColumnIndex]);

    const metricHISASPColumnIndex = fetchData.dict.columns.indexOf('metric_HIS-ASP');
    const metricHISASPColumn = fetchData.dict.data.map(row => row[metricHISASPColumnIndex]);

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
    
    const distances = fetchData.distances[proteinSelector.value][ligandSelector.value]

    console.log(distances)
    
    // Combina columnNames con las distancias dinámicas
    const allColumnNames = [...columnNames, ...distances];
    
    console.log(allColumnNames)

    // Mapea los nombres de las columnas a las definiciones de columnas de ag-Grid
    const columnDefs = allColumnNames.map(column => {
        return { field: column };
    });
    
    const rowData = proteinColumn.map((protein, index) => {
        return {
            "Protein": protein,
            "Ligand": ligandColumn[index],
            "Epoch": epochColumn[index],
            "Trajectory": trajectoryColumn[index],
            "Accepted Pele Steps": acceptedPeleStepsColumn[index],
            "Step": stepColumn[index],
            "Total Energy": totalEnergyColumn[index],
            "Binding Energy": bindingEnergyColumn[index],
            "Ligand SASA": ligandSASAColumn[index],
            "metric_SER-L": metricSERLColumn[index],
            "metric_SER-HIS": metricSERHISColumn[index],
            "metric_HIS-ASP": metricHISASPColumn[index]
        };
    });

    distances.forEach((distance, distanceIndex) => {
        const distanceColumnIndex = 12 + distanceIndex;
        const distanceColumn = fetchData.dict.data.map(row => row[distanceColumnIndex]);
    
        rowData.forEach((row, rowIndex) => {
            row[distance] = distanceColumn[rowIndex];
        });
    });

    const gridOptions = {
        rowData: rowData,
        columnDefs: columnDefs,
        suppressSizeToFit: true
    };
    
    const gridElement = document.querySelector('#myGrid');
    const gridApi = agGrid.createGrid(gridElement, gridOptions);
}

const showPlot = () => {
    let distanceValue = distanceSelector.value
    const distanceColumnIndex = fetchData.dict.columns.indexOf(distanceValue);
    distanceColumn = fetchData.dict.data.map(row => row[distanceColumnIndex]);

    let colorSelectorIndex = fetchData.dict.columns.indexOf(colorSelector.value);
    let colorSelectorColumn = fetchData.dict.data.map(row => row[colorSelectorIndex]);

    showGrid()

    if(colorSelector.value == 'None'){
        var pelePlot = document.getElementById('pelePlot'),
        data = [{
            x: distanceColumn, 
            y: dataColumn, 
            type: 'scatter',
            mode: 'markers', 
            marker: {
                size: 6,
                color: 'rgba(45, 85, 255, 0.7)',    
            }
                }],
        layout = {
            title: 'Scatter plot',
            hovermode: 'closest',
            width: 500,
            height: 500,
            xaxis: {
                title: distanceValue
            },
            yaxis: {
                title: 'Binding energy'
            },
            marker: { size: 3 }
        };
    
        Plotly.newPlot('pelePlot', data, layout);
    }else{
        let trace = {
            x: distanceColumn,
            y: dataColumn,
            mode: 'markers',
            type: 'scatter',
            marker: {
                size: 4,
                color: colorSelectorColumn,
                colorscale: 'Viridis',
                colorbar: {
                    title: colorSelector.value
                }
            },
        };
        
        let layout = {
            title: 'Scatter plot',
            width: 600,
            height: 500,
            xaxis: {
                title: distanceValue
            },
            yaxis: {
                title: 'Binding energy'
            },
            
        };
        Plotly.newPlot('pelePlot', [trace], layout);
    }
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
            .then(response => {
                response.json()
            })
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

                    showPlot()

                    pelePlot.on('plotly_click', function(){
                        alert('You clicked this Plotly chart!');
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

                    showPlot()
                }else{
                    console.log('Error: ' + parsedData.msg)
                }
            })
            .catch(error => console.error('Error en la solicitud:', error));
    }


    showPlot()

    pelePlot.on('plotly_click', function(){
        alert('You clicked this Plotly chart!');
    });
})
