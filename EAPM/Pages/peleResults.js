const peleSimulationFolder = document.querySelector('#peleSimulationFolder')
const peleSimulationFolderBtn = document.querySelector('#peleSimulationFolderBtn')
const peleOutputFolder = document.querySelector('#peleOutputFolder')

peleSimulationFolderBtn.addEventListener("click", (e) => {
    e.preventDefault()

    const href = window.location.href

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
            console.log(parsedData);
        })
        .catch(error => console.error('Error en la solicitud:', error));
})

