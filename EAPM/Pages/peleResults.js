const peleSimulationFolder = document.querySelector('#peleSimulationFolder')
const peleSimulationFolderBtn = document.querySelector('#peleSimulationFolderBtn')
const peleOutputFolder = document.querySelector('#peleOutputFolder')
const peleForm = document.querySelector('#peleForm')
const plotForm = document.querySelector('#plotForm')
const proteinSelector = document.querySelector('#proteinSelector')
const ligandSelector = document.querySelector('#ligandSelector')
const distanceSelector = document.querySelector('#distanceSelector')
const colorSelector = document.querySelector('#colorSelector')

let fetchData
let dataColumn

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

                showPlot()
            }else{
                console.log('Error: ' + parsedData.msg)
            }
        })
        .catch(error => console.error('Error en la solicitud:', error));
})

const showPlot = () => {
    let x = fetchData.dict.columns[6];
    let y = fetchData.dict.columns[7];

    let distanceValue = distanceSelector.value
    const distanceColumnIndex = fetchData.dict.columns.indexOf(distanceValue);
    distanceColumn = fetchData.dict.data.map(row => row[distanceColumnIndex]);

    let colorSelectorIndex = fetchData.dict.columns.indexOf(colorSelector.value);
    let colorSelectorColumn = fetchData.dict.data.map(row => row[colorSelectorIndex]);

    if(colorSelector.value == 'None'){
        let trace = {
            x: distanceColumn,
            y: dataColumn,
            mode: 'markers',
            type: 'scatter',
            marker: {
                size: 6,
                color: 'rgba(45, 85, 255, 0.7)',    
            }
        };
    
        let layout = {
            title: 'Scatter plot',
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
    
        Plotly.newPlot('testPlot', [trace], layout);
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
        Plotly.newPlot('testPlot', [trace], layout);
    }
}

plotForm.addEventListener("change", (e) => {

    if(e.target.id === "proteinSelector"){

        const data = {
            peleSimulationFolder: peleSimulationFolder.value,
            peleOutputFolder: peleOutputFolder.value,
            desiredProtein: proteinSelector.value,
            desiredLigand: ligandSelector.value
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
    
                    // const proteins = parsedData.distances
                    // const proteinValue = Object.keys(proteins)[0]
                    // const ligands = parsedData.distances[proteinValue]
                    // const distances = parsedData.distances[proteinValue][Object.keys(ligands)[0]]
                    
                    ligandSelector.innerHTML = ''
                    distanceSelector.innerHTML = ''

                    const proteinValue = e.target.value
                    const ligands = parsedData.distances[proteinValue]
                    const distances = parsedData.distances[proteinValue][ligandSelector.value]

                    addDictOption(ligandSelector, ligands)

                    addOption(distanceSelector, distances)
                    
                    showPlot()
                }else{
                    console.log('Error: ' + parsedData.msg)
                }
            })
            .catch(error => console.error('Error en la solicitud:', error));
    }


    showPlot()
})



