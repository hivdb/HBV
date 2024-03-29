document.getElementById('csvFileInput').addEventListener('change', loadFile);

function loadFile(event) {
    const input = event.target;
    if ('files' in input && input.files.length > 0) {
        readFileContent(input.files[0]).then(content => {
            const data = parseCSV(content);
            displayTable(data);
            document.getElementById('csvFileInput').style.display = 'none';
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
        var overall_cons = rows[0].overall_cons;

        if (index % num_per_row === 0) {
            var cell = document.createElement('td');
            var posElement = document.createElement('div')
            posElement.textContent = ''
            cell.appendChild(posElement);
            table.lastChild.appendChild(cell);
        }

        var cell = document.createElement('td');
        cell.className = 'pos';
        var posElement = document.createElement('div');
        index += 1
        if ((index % num_per_row) === 0 || (index % 5 == 0) || (index === 1)) {
            posElement.textContent = pos;
        } else {
            posElement.textContent = ' '
        }
        var consElement = document.createElement('div');
        consElement.textContent = overall_cons;
        cell.appendChild(posElement);
        // cell.appendChild(consElement);
        table.lastChild.appendChild(cell);
    });

    data = data.map((pos) => pos[1]).flat()

    data = groupByGenotype(data)

    const sortedKeys = Object.keys(data).sort();
    var sortedGroup = sortedKeys.map(key => [key, data[key]]);

    var lastElement = sortedGroup.pop();
    sortedGroup.unshift(lastElement);

    sortedGroup.forEach(([genotype, row]) => {
        displayGenotype(genotype, row, table)
    })

}

function displayGenotype(genotype, data, table) {

    groupData = groupByPos(data);

    table.appendChild(document.createElement('tr'));

    var cell = document.createElement('td');
    cell.className = 'genotype';

    var posElement = document.createElement('div')
    posElement.textContent = genotype + " (" + data[0].total + ")";
    cell.appendChild(posElement);
    table.lastChild.appendChild(cell)

    Object.keys(groupData).sort((a, b) => a - b).forEach((pos, index) => {
        var cell = document.createElement('td');
        cell.className = 'mutationGroup';

        groupData[pos].sort((a, b) => b.pcnt - a.pcnt).forEach(mut => {
            const mutElement = document.createElement('div');
            mutElement.textContent = mut.mut;
            const pcntElement = document.createElement('sup');
            pcntElement.textContent = mut.pcnt;
            mutElement.className = 'mutations';
            if (mut.pcnt >= 90 && mut.mut == mut.overall_cons) {
                mutElement.className += ' pcnt100';
            } else if (mut.pcnt >= 50) {
                mutElement.className += ' pcnt50';
            } else if (mut.pcnt >= 10) {
                mutElement.className += ' pcnt10';
            } else if (mut.pcnt >= 1) {
                mutElement.className += ' pcnt1';
            } else {
                mutElement.className += ' pcnt01';
            }

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
