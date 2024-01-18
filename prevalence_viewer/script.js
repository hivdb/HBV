document.getElementById('csvFileInput').addEventListener('change', loadFile);

function loadFile(event) {
    const input = event.target;
    if ('files' in input && input.files.length > 0) {
        readFileContent(input.files[0]).then(content => {
            const data = parseCSV(content);
            displayTable(data);
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
    const lines = text.split(/\r\n|\n/).filter(string => string !== '');
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


function displayTable(data) {
    const tableContainer = document.getElementById('tableContainer');
    tableContainer.innerHTML = ''; // Clear previous table

    const table = document.createElement('table');

    const posGroups = groupByPos(data);
    const sortedKeys = Object.keys(posGroups).sort((a, b) => a - b);

    num_per_row = 50

    const sortedPosGroups = sortedKeys.map(key => [key, posGroups[key]]);
    const batches = splitArrayIntoChunks(sortedPosGroups, num_per_row);

    batches.forEach((row) => {
        displayRow(row, table)
    })

    tableContainer.appendChild(table);
}

function displayRow(data, table) {

    table.appendChild(document.createElement('tr'));

    data.forEach((pair, index) => {
        var pos = pair[0];
        var rows = pair[1];

        if (index % num_per_row === 0) {
            var cell = document.createElement('td');
            var posElement = document.createElement('div')
            posElement.className = 'pos'
            posElement.textContent = 'Genotype'
            cell.appendChild(posElement);
            table.lastChild.appendChild(cell);
        }

        var cell = document.createElement('td');
        var posElement = document.createElement('div');
        posElement.className = 'pos';
        posElement.textContent = pos;
        cell.appendChild(posElement);
        table.lastChild.appendChild(cell);
    });

    data = data.map((pos) => pos[1]).flat()

    data = groupByGenotype(data)

    const sortedKeys = Object.keys(data).sort();
    const sortedGroup = sortedKeys.map(key => [key, data[key]]);

    sortedGroup.forEach(([genotype, row]) => {
        displayGenotype(genotype, row, table)
    })

}

function displayGenotype(genotype, data, table) {

    data = groupByPos(data);

    table.appendChild(document.createElement('tr'));

    var cell = document.createElement('td');
    var posElement = document.createElement('div')
    posElement.textContent = genotype;
    cell.appendChild(posElement);
    table.lastChild.appendChild(cell)

    Object.keys(data).sort((a, b) => a - b).forEach((pos, index) => {
        var cell = document.createElement('td');

        data[pos].sort((a, b) => b.pcnt - a.pcnt).forEach(mut => {
            const mutElement = document.createElement('div');
            mutElement.textContent = mut.mut;
            const pcntElement = document.createElement('sup');
            pcntElement.textContent = mut.pcnt;
            mutElement.className = 'mutations';
            mutElement.appendChild(pcntElement);
            cell.appendChild(mutElement);
        });


        // cell.children[1].className = 'consensus'

        table.lastChild.appendChild(cell);
    });

}

function groupByPos(data) {
    return data.reduce((groups, item) => {
        let pos = parseInt(item.pos);
        if (!groups[pos]) {
            groups[pos] = [];
        }
        groups[pos].push(item);
        return groups;
    }, {});
}


function groupByGenotype(data) {
    return data.reduce((groups, item) => {
        if (!groups[item.genotype]) {
            groups[item.genotype] = [];
        }
        groups[item.genotype].push(item);
        return groups;
    }, {});
}


function splitArrayIntoChunks(array, chunkSize) {
    let result = [];
    for (let i = 0; i < array.length; i += chunkSize) {
        let chunk = array.slice(i, i + chunkSize);
        result.push(chunk);
    }
    return result;
}
