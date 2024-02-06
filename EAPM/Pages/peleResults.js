const peleSimulationFolder = document.querySelector('#peleSimulationFolder')
const peleSimulationFolderBtn = document.querySelector('#peleSimulationFolderBtn')
const peleOutputFolder = document.querySelector('#peleOutputFolder')
const peleForm = document.querySelector('#peleForm')
const plotForm = document.querySelector('#plotForm')
const plotButton = document.querySelector('#plotButton')
const proteinSelector = document.querySelector('#proteinSelector')
const ligandSelector = document.querySelector('#ligandSelector')
const distanceSelector = document.querySelector('#distanceSelector')
const colorSelector = document.querySelector('#colorSelector')

let fetchData

const addOption = (selector, array) => {
    array.forEach(i => {
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

    let result = {}

    const data = {
        peleSimulationFolder: peleSimulationFolder.value,
        peleOutputFolder: peleOutputFolder.value
    }

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

                const proteins = parsedData.proteins
                const ligands = parsedData.ligands
                const distances = parsedData.distances

                addOption(proteinSelector, proteins)
                addOption(ligandSelector, ligands)
                addOption(distanceSelector, distances)

                peleForm.style.display = 'none'

                plotForm.style.display = 'block'
            }
        })
        .catch(error => console.error('Error en la solicitud:', error));
})

plotButton.addEventListener("click", (e) => {
    e.preventDefault()

    let x = fetchData.dict.columns[6];
    let y = fetchData.dict.columns[7];

    let distanceValue = distanceSelector.value

    const totalEnergy = fetchData.dict.columns.indexOf("Binding Energy");
    let colorSelectorIndex = fetchData.dict.columns.indexOf(colorSelector.value);
    const distanceColumn = fetchData.dict.columns.indexOf(distanceValue);

    const columnData = fetchData.dict.data.map(row => row[totalEnergy]);
    const columnDistance = fetchData.dict.data.map(row => row[distanceColumn]);
    let colorSelectorColumn = fetchData.dict.data.map(row => row[colorSelectorIndex]);


    if(colorSelector.value == 'None'){
        let trace = {
            x: columnData,
            y: columnDistance,
            mode: 'markers',
            type: 'scatter',
            marker: {
                size: 10,
                color: 'rgba(50, 171, 96, 0.7)',
            }
        };
    
        let layout = {
            title: 'Scatter plot',
            xaxis: {
                title: distanceValue
            },
            yaxis: {
                title: 'Binding energy'
            }
        };
    
        Plotly.newPlot('testPlot', [trace], layout);
    }else{
        let trace = {
            x: columnData,
            y: columnDistance,
            mode: 'markers',
            type: 'scatter',
            marker: {
                size: 10,
                color: colorSelectorColumn,
                colorscale: 'Viridis',
                colorbar: {
                    title: colorSelector.value
                }
            }
        };
        
        let layout = {
            title: 'Scatter plot',
            xaxis: {
                title: distanceValue
            },
            yaxis: {
                title: 'Binding energy'
            }
        };
        
        Plotly.newPlot('testPlot', [trace], layout);
    }
})



