function importModule (dir) {
    let mod = document.createElement("script");
    mod.src = dir;
    mod.type = "text/javascript";
    document.head.appendChild(mod);
}

importModule("/calc.js")

//Macros
const maxRows = 512;
const PI = Math.PI;
const TWOPI = 2 * PI;
const E = Math.E;
const sampleRate = 44100;
const roundCutoff = 7;

function round(num) {
    let roundCutoffNum = Math.pow(10, roundCutoff) * 1.0;
    return Math.round(roundCutoffNum * num) * 1.0 / roundCutoffNum;
}


function userFunc (func, table) {

    let direction = (func.includes('i')) ? -1 : 1;

    let arg = collectSamples(table, direction);

    // console.log(arg);
    let oldSize = arg.real.length;
    let nextPow2 = Math.pow(2, Math.ceil(Math.log2(oldSize)));
    // Pad with 0s
    while (arg.real.length < nextPow2) {
        arg.real.push(0);
        arg.imaginary.push(0);
    }

    // If we don't already have imaginaries
    if (direction == 1)
        for (let i = 0; i < arg.real.length; i++)
            arg.imaginary[i] = 0; 

    let result;

    if (func.includes('fft'))
        result = fft(arg, direction);
    else if (func.includes('dft'))
        result = dft(arg, direction);


    printResults(result)
    
}

function printResults(result) {
    let resultTable = document.querySelector("#resultsTable");
    setRowNum(resultTable, 0)
    for (let i = 0; i < result.real.length; i++) {
        let r = result.real[i];
        let row = addRow(resultTable);
        let c1 = row.cells[1];
        let c2 = row.insertCell(-1);
        let c3 = row.insertCell(-1);
        c2.innerHTML = c1.innerHTML;
        c3.innerHTML = c1.innerHTML;
        let cells = [c1, c2, c3]
        let values = [];
        let j = 0;
        cells.forEach(cell=>{values.push(cell.getElementsByTagName("input")[0])});
        values[0].value = round(result.real[i]);
        values[1].value = round(result.imaginary[i]);
        values[2].value = round(Math.sqrt(result.real[i]*result.real[i] + result.imaginary[i]*result.imaginary[i]));
    }
}

function addRow(table) {

    if (table.rows.length > maxRows)
        return;

    let row = table.insertRow(-1);
    let c1 = row.insertCell(0);
    c1.innerHTML = (table.rows.length - 1);
    c1.style.textAlign="center";
    row.insertCell(1).innerHTML = "<input type='text' class='txt' placeholder='Type to enter value'/>";
    return row;
}

function deleteRow(table) {

    if (table.rows.length > 1)
        table.deleteRow(-1);   
}

function setRowNum(table, n) {
    n++;

    if (n > 0 && n < table.rows.length)
        while (n < table.rows.length)
            deleteRow(table);
    else if (n <= maxRows + 1 && n > table.rows.length)
        while (n > table.rows.length)
            addRow(table);
}

function collectSamples(table, direction) {
    let i = 0;
    let samples = {};
    samples.real = [];
    samples.imaginary = [];
    Array.from(table.rows).forEach(row=>{
        if (i++ != 0) {
            samples.real[i-2] = parseFloat(row.cells[1].getElementsByTagName("input")[0].value);
            if (direction == -1)
                samples.imaginary[i-2] = parseFloat(row.cells[2].getElementsByTagName("input")[0].value);
            else
                samples.imaginary.push(0)
        }
    });
    return samples;
}

function fillRandom(table) {
    let i = 0;
    Array.from(table.rows).forEach(row=>{
        if (i++ != 0) {
            row.cells[1].getElementsByTagName("input")[0].value = Math.random();
        }
    });
}

function loadD1(table) {

    let samples = [ 0, 1, 2, 3, 4, 5, 6, 7 ];
    setRowNum(table, 0)
    setRowNum(table, 8)

    let i = 0;
    Array.from(table.rows).forEach(row=>{
        if (i++ != 0) {
            row.cells[1].getElementsByTagName("input")[0].value = samples[i-2];
        }
    });
}

function loadD2(table) {

    let samples = [];

    for (let i = 0; i < maxRows; i++)
        samples[i] = Math.sin(TWOPI * 244 * i / maxRows )

    setRowNum(table, 0)
    setRowNum(table, maxRows)

    let i = 0;
    Array.from(table.rows).forEach(row=>{
        if (i++ != 0) {
            row.cells[1].getElementsByTagName("input")[0].value = samples[i-2];
        }
    });
}

document.addEventListener("load", ()=>{init()});