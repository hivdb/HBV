document.getElementById('csvFileInput').addEventListener('change', loadFile);

function loadFile(event) {
    const input = event.target;
    if ('files' in input && input.files.length > 0) {
        readFileContent(input.files[0]).then(content => {
            const data = parseCSV(content);
            populateDropdown(data);
        }).catch(error => console.log(error));
    }
}

function readFileContent(file) {
    const reader = new FileReader();
    return new Promise((resolve, reject) => {
        reader.onload = event => resolve(event.target.result);
        reader.onerror = error => reject(error);
        reader.readAsText(file);
    });
}

function parseCSV(text) {
    const lines = text.split(/\r\n|\n/);
    const headers = lines[0].split(',');
    const data = lines.slice(1).map(line => {
        const values = line.split(',');
        return headers.reduce((object, header, index) => {
            object[header] = values[index];
            return object;
        }, {});
    });
    return data;
}

function populateDropdown(data) {
    const genotypeSet = new Set(data.map(item => item.genotype));
    const dropdown = document.getElementById('genotypeDropdown');
    genotypeSet.forEach(genotype => {
        const option = document.createElement('option');
        option.value = genotype;
        option.textContent = genotype;
        dropdown.appendChild(option);
    });

    dropdown.addEventListener('change', () => {
        const selectedGenotype = dropdown.value;
        const filteredData = data.filter(item => item.genotype === selectedGenotype);
        displayTable(filteredData);
    });
}

function displayTable(data) {
    const tableContainer = document.getElementById('tableContainer');
    tableContainer.innerHTML = ''; // Clear previous table

    const posGroups = groupByPos(data);
    const table = document.createElement('table');

    Object.keys(posGroups).sort((a, b) => a - b).forEach((pos, index) => {
        if (index % 50 === 0) { // Start a new row every 50 positions
            table.appendChild(document.createElement('tr'));
        }

        const cell = document.createElement('td');
        posGroups[pos].sort((a, b) => b.pcnt - a.pcnt).forEach(mut => {
            const mutElement = document.createElement('div');
            mutElement.textContent = mut.mut + ' ';
            const pcntElement = document.createElement('sup');
            pcntElement.textContent = mut.pcnt;
            mutElement.appendChild(pcntElement);
            cell.appendChild(mutElement);
        });

        table.lastChild.appendChild(cell);
    });

    tableContainer.appendChild(table);
}

function groupByPos(data) {
    return data.reduce((groups, item) => {
        if (!groups[item.pos]) {
            groups[item.pos] = [];
        }
        groups[item.pos].push(item);
        return groups;
    }, {});
}
