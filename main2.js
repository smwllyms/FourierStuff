function importModule (dir) {
    let mod = document.createElement("script");
    mod.src = dir;
    mod.type = "text/javascript";
    document.head.appendChild(mod);
}

let barsElem;
let detectedPitchElem = document.getElementById("detectedPitch");
let freqBars = [];

const freqOfMiddleC = 256;
// credit: https://newt.phys.unsw.edu.au/jw/notes.html#:~:text=at%20the%20top.-,To%20convert%20from%20any%20frequency%20to%20pitch%20(i.e.%20to%20the,converter%20written%20by%20Andrew%20Botros.&text=no%20%3D%20log2(f2%2Ff1).&text=fn%20%3D%202n%2F12*440%20Hz.
const getPitch = function(hertz) {
    // Let's be a little sharper because with this alg below C is B
    // let offset = 5;
    let pitchOffsetSemitones = 12*Math.log2(hertz/freqOfMiddleC)
    // corrected, my switch starts at A which is 2 steps below
    let corrected = Math.floor(pitchOffsetSemitones);
    let octaveNum = Math.floor(Math.log2(hertz/freqOfMiddleC)) + 5;
    if (pitchOffsetSemitones >= 0)
        corrected = corrected % 12;
    else {
        corrected = 12 - (-corrected % 12);
    }
    switch (corrected) {
        case 9:
            return "A" + octaveNum;
        case 10:
            return "A#/Bb" + octaveNum;
        case 11:
            return "B" + octaveNum;
        case 0:
            return "C" + octaveNum;
        case 1:
            return "C#/Db" + octaveNum;
        case 2:
            return "D" + octaveNum;
        case 3:
            return "D#/Eb" + octaveNum;
        case 4:
            return "E" + octaveNum;
        case 5:
            return "F" + octaveNum;
        case 6:
            return "F#/Gb" + octaveNum;
        case 7:
            return "G" + octaveNum;
        case 8:
            return "G#/Ab" + octaveNum;
        default:
            return "A" + octaveNum;
    }
}

const drawBars = function(data) {
    let height = parseInt(barsElem.offsetHeight) - 5;
    let effectiveHeight;
    // More accurate prediction
    let cutoffFreq = 4000;
    let effectiveLength = cutoffFreq / 44100 * data.length;
    // let max = Math.max(...data);
    let max = -999, average = 0, totalBias = 0, totalNoBias = 0;
    let maxI = 0;
    for (let i = 0; i < effectiveLength; i++) {
        if (data[i] > max) {
            max = data[i];
            maxI = i;
        }
        average += data[i];
    }
    average /= effectiveLength;
    do {
        let i = maxI;
        let dif = 1;
        while (data[i] >= average) {
            dif = Math.abs(data[i] - data[i+1])/data[i+1];
            totalNoBias += data[i] * dif;
            totalBias += i * data[i] * dif;
            i--;
        }
        i = maxI + 1;
        while (data[i] >= average) {
            dif = Math.abs(data[i] - data[i-1])/data[i-1];
            totalNoBias += data[i] * dif;
            totalBias += i * data[i] * dif;
            i++;
        }
    } while (false);

    let expectedBin = totalBias / totalNoBias;
    console.log((expectedBin/ (data.length * 2) * 44100) + ", " + (maxI/ (data.length * 2) * 44100))

    // Color rectangles
    let bin = 0, val;
    while (bin < effectiveLength) {
        val = data[bin]
        effectiveHeight = (val / max * height);
        if (effectiveHeight < 1)
            effectiveHeight = 1;
        else if (effectiveHeight > height)
            effectiveHeight = height;
        let binPos = bin / data.length * 8 * 255;
        freqBars[bin].style.backgroundColor = "rgba(" + binPos + ", 128, " + (255 - binPos) + ", 1.0)";
        freqBars[bin++].style.height = effectiveHeight + "px";
    }
    // My special twiddle factor 1.087 from desmos graphing expected vs actual
    let guess = Math.round(expectedBin / (data.length * 2) * 44100 * 1.087);
    let pitch = getPitch(guess);
    detectedPitchElem.innerHTML = "Detected Pitch: "+ guess + " Hz (" + pitch + ")";   
}


let start = document.getElementById("start");
let form = document.getElementsByTagName("form")[0];
let freqSlider;
let actualOscPitch = document.getElementById("actualOscPitch");
let mode = 1; // 1: sine, 2: mic
let running = false;

start.onclick = function() { 
    if (!running) {
        form.dispatchEvent(new Event("submit")); 
        running = true; 
    }
    else location.reload();
}

form.onsubmit = async function(e) {
    e.preventDefault();
    const formData = new FormData(e.target);
    const formProps = Object.fromEntries(formData);
    
    if (formProps.mode)
        mode = formProps.mode;

    // if (mode == 1)
    
    const audioContext = new AudioContext();
    await audioContext.audioWorklet.addModule('worker.js')
    const workerNode = new AudioWorkletNode(audioContext, 'worker-node', {processorOptions:{size:1024}});
    // ON messgae paint bars
    workerNode.port.onmessage = (e) => drawBars(e.data);
    let streamSrc;

    if (mode == 3) {
        await navigator.mediaDevices.getUserMedia({ audio: true, video: false })
        .then(stream=>streamSrc = audioContext.createMediaStreamSource(stream));
    }
    else {
        let osc = audioContext.createOscillator();
        let scale = 16, freq=freqOfMiddleC;
        osc.frequency.value = freq;
        // Saw or sine
        switch (parseInt(mode)) {
            case 1:
                osc.type = "sine";
                break;
            case 2:
                osc.type = "sawtooth";
                break;
            default:
                osc.type = "sine";
                break;
        }
        osc.start();
        streamSrc = audioContext.createGain();
        streamSrc.gain.value = 0.1;
        osc.connect(streamSrc);
        freqSlider.addEventListener("input", (e)=>{
            osc.frequency.value = freq * Math.pow(2, 2*(e.target.value - 50) / 100); 
            actualOscPitch.innerText = getPitch(Math.round( 100 * osc.frequency.value) / 100) + " (" + Math.round( 100 * osc.frequency.value) / 100 + " Hz)";
        });
    }
    streamSrc.connect(workerNode);
    workerNode.connect(audioContext.destination)
}

// ON document load
barsElem = document.getElementById("displayBars");
freqSlider = document.getElementById("freqSlider");
for (let i = 0; i < 1024; i++) {
    let bar = document.createElement("div");
    bar.classList.add("bar")
    barsElem.appendChild(bar);

    freqBars[i] = bar;
}