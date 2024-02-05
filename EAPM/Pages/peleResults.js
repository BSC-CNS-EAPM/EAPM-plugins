const peleSimulationFolder = document.querySelector('#peleSimulationFolder')
const peleSimulationFolderBtn = document.querySelector('#peleSimulationFolderBtn')
const peleOutputFolder = document.querySelector('#peleOutputFolder')
const peleForm = document.querySelector('#peleForm')
const plotForm = document.querySelector('#plotForm')
const proteinSelector = document.querySelector('#proteinSelector')
const ligandSelector = document.querySelector('#ligandSelector')
const distanceSelector = document.querySelector('#distanceSelector')

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
                console.log(parsedData)

                // TODO Optimizar el código repetitivo

                const proteins = parsedData.proteins
                const ligands = parsedData.ligands
                const distances = parsedData.distances

                proteins.forEach(i => {
                    protein = document.createElement('option')
                    protein.id = i
                    protein.value = i
                    protein.innerHTML = i
                    proteinSelector.append(protein)
                });

                ligands.forEach(i => {
                    ligand = document.createElement('option')
                    ligand.id = i
                    ligand.value = i
                    ligand.innerHTML = i
                    ligandSelector.append(ligand)
                });

                distances.forEach(i => {
                    distance = document.createElement('option')
                    distance.id = i
                    distance.value = i
                    distance.innerHTML = i
                    distanceSelector.append(distance)
                });

                peleForm.style.display = 'none'

                plotForm.style.display = 'block'

                var data = [{
                    values: [19, 26, 55],
                    labels: ['Test1', 'Test2-', 'Test3'],
                    type: 'pie'
                  }];
                  
                  var layout = {
                    height: 400,
                    width: 500
                  };
                  
                  Plotly.newPlot('testPlot', data, layout);
            }
        })
        .catch(error => console.error('Error en la solicitud:', error));
})

