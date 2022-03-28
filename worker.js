class WorkerNode extends AudioWorkletProcessor {

    constructor(props) {
        super()
        this.idx = 0;
        this.size = props.processorOptions.size;
        console.log(this.size)
        this.buffer = new Float32Array(this.size);
        this.sampleRate = props.processorOptions.sampleRate;

        this.window = [];
        for (let i = 0; i < this.size; i++) {
            // Hahn window should improve accuracy of FT
            this.window[i] = 0.5 - 0.5*Math.cos((2*Math.PI*i)/(this.size-1));
        }
    }

    fft(inReal, inImaginary, inv) {
        let tr, ti;
        let len = this.size;
        let arr = [];
        let M_PI = Math.PI;
        let direction = inv;

        // 1. Swap bits
        let logSize = Math.log2(len);
        for (let i = 0; i < len; i++) {
            let reversed = 0;
            for (let j = 0; j < logSize; j++) {
                reversed |= !!((1 << j) & i) << (logSize - j - 1);
            }
            // console.log(i, reversed)
            // Just so we only change once
            if (reversed >= i) {
    
                let tmp = inReal[i];
                inReal[i] = inReal[reversed];
                inReal[reversed] = tmp;
    
                tmp = inImaginary[i];
                inImaginary[i] = inImaginary[reversed];
                inImaginary[reversed] = tmp;
            }
        }

        // 1.5 assign 2 arrays to 1
        for (let i = 0; i < len; i++) {
            arr[i] = inReal[i] * this.window[i];
            arr[i+len] = inImaginary[i] * this.window[i];
        }
        let size = len;
        let TWOPI = 2*M_PI;

        // 2. fft
        for (let i = 1; i <= Math.log2(len); i++) {
            let le = Math.pow(2, i);
            let le2 = le / 2;
            let ur = 1;
            let ui = 0;
            let sr = Math.cos(M_PI / le2);
            let si = direction * -Math.sin(M_PI / le2);
            for (let j = 0; j < le2; j++) {
                for (let k = j; k < len; k += le) {
                    let kp = k + le2;
                    tr = arr[kp] * ur - arr[kp + len] * ui;
                    ti = arr[kp] * ui + arr[kp + len] * ur;
                    arr[kp] = arr[k] - tr;
                    arr[kp + len] = arr[k + len] - ti;
                    arr[k] += tr;
                    arr[k + len] += ti;
                }
                tr = ur;
                ur = tr * sr - ui * si;
                ui = tr * si + ui * sr;
            }
        }
        for (let i = 0; i < size; i++) {
            inReal[i] = arr[i];
            inImaginary[i] = arr[i+len];
        }
        if (inv == -1) {
            for (let i = 0; i < len; i++) {
                inReal[i] /= len;
            }
        }
    }

    processBuffer() {
        /* Manipulate the sample */
        let real = this.buffer
        let imaginary = new Float32Array(real.length);
        this.fft(real, imaginary, 1);
        // We have freq data

        // Get magnitude
        let mag = new Float32Array(real.length/2);
        for (let i = 0; i < real.length/2; i++)
            mag[i] = Math.sqrt(real[i]*real[i]+imaginary[i]*imaginary[i])
        // 1. display data
        this.port.postMessage(mag);

        // Create new time data
        // this.fft(real, imaginary, -1);
    }

    process(inputList, outputList, parameters) {
        const sourceLimit = Math.min(inputList.length, outputList.length);
      
        for (let inputNum = 0; inputNum < sourceLimit; inputNum++) {
          let input = inputList[inputNum];
          let output = outputList[inputNum];
          let channelCount = Math.min(input.length, output.length);
      
          for (let channelNum = 0; channelNum < channelCount; channelNum++) {
                let sampleCount = input[channelNum].length;
      
                output[channelNum].set(input[channelNum]);
                this.buffer.set(input[channelNum], this.idx);

                this.idx += sampleCount;
                if (this.idx == this.size) {
                    this.processBuffer();
                    this.idx = 0;
                }
          }
        };
      
        return true;
      }
}

registerProcessor('worker-node', WorkerNode)