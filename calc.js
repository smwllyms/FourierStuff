function dft(samples, direction) {
    if (samples.real.length < 2) 
        return;

    let inReal = samples.real;
    let inImaginary = samples.imaginary;
    let outReal = [];
    let outImaginary = [];
    let size = inReal.length;
    for (let i = 0; i < size; i++) {
        outReal[i] = 0;
        outImaginary[i] = 0;
    }

    let time = new Date().getTime();
    // For each value of waveform, consider real and imaginary components
    // e^(aix) = cos(ax)+isin(ax)
    // in our case, real is cos(ax) and imaginary is isin(ax) (really just holds sin(ax))
    // DFT = sum k=1 -> size [ g(tk) * e^(-2pi*f*tk) ], where in our case tk = k
    // * we start on 0 instead of 1 and f = bin / (sample size)
    // * gtk is the value of the sample at time tk = gtkr, gtki is imaginary used for inverse
    // * imaginary is negative (convention is clockwise)
    let gtkr, gtki, a, f;
    for (let bin = 0; bin < size; bin++) {
        // calculate frequency
        f = bin / size;
        // Coefficient
        a = -TWOPI * f * direction;
        // For each sample (we are summing all samples)
        for (let k = 0; k < size; k++) {
            gtkr = inReal[k]; // real sample
            gtki = inImaginary[k]; // imaginary sample (not used for pure time input)
            // out real bin is sum of real (cos) part (imaginary is 0 for pure time input)
            outReal[bin] += gtkr * Math.cos(a * k) + gtki * Math.sin(a * k); 
            // out imaginary bin is sum of imaginary (sin) (imaginary is 0 for pure time input)
            outImaginary[bin] += -(gtkr * Math.sin(a * k) + gtki * Math.cos(a * k)); // negative from clockwise
        }
    }

    console.log("DFT Completed in " + (new Date().getTime() - time) + " milliseconds!")

    // If inverse we have to divide by the size of input (our reference)
    if (direction == -1) {
        for (let i = 0; i < size; i++) {
            outReal[i] /= size;
            outImaginary[i] /= size;
        }
    }

    // Now we apply to out table
    let result = {};
    result.real = outReal;
    result.imaginary = outImaginary;
    return result;
}
let mode = -1;
function fft(samples, direction) {
    if (samples.real.length < 2) 
        return;

    let inReal = samples.real;
    let inImaginary = samples.imaginary;
    let outReal = [];
    let outImaginary = [];
    let size = inReal.length;
    for (let i = 0; i < size; i++) {
        outReal[i] = 0;
        outImaginary[i] = 0;
    }

    let time = window.performance.now();

    // First we must flip the bits (think repeatedly take odds and evens)
    // console.log(inReal)
    let logSize = Math.log2(size);
    for (let i = 0; i < size; i++) {
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
    // console.log(inReal)

    mode = (mode + 1) % 4;
    if (mode == 0) {
        // Radix 2
        // Log2 size stages (see diagram for how many butterfly stages there are)
        for (let stage = 1; stage <= logSize; stage++) {
            // Separation between butterfly chains = total length of each butterfly chain
            let separation = Math.pow(2, stage);
            // The width of each butterfly chain = distance between lo and hi in DFT
            let width = separation / 2;
            // Twiddle factors for identifying where on Euler's identity circle we are 
            // Rotation in complex coordinates. 
            //
            // We will rotate/increment twiddle factors by 1 / separation counter-clockwise
            // each stage by using the formula rot (a rads) = base # * e^( i * a )
            //
            // So in programming terms:
            // rot (val, rads): 
            //      result.real = val.real * [ cos( a ) ] + val.imaginary * [ sin (a) ]
            //      result.imaginary = -val.real * [ sin( a ) ] + val.imaginary * [ cos (a) ]
            //
            // Start on positive, real axis (1+0j), (aka 0deg)
            let twiddleReal = 1;
            let twiddleImaginary = 0;
            // For each offset of butterfly chain
            for (let i = 0; i < width; i++) {
                // ... in each butterfly chain
                for (let currentButterfly = i; currentButterfly < size; currentButterfly += separation) {

                    // Get hi and lo ptrs for DFT
                    let hi = currentButterfly;
                    let lo = hi + width;

                    // Algorithm states to take DFT of lo and add it to hi, and negate
                    // the value at lo and add the hi

                    // 1. Get value of lo (rel and imaginary)
                    // Recall DFT = sum_k=0^N ( g(t_k) * e^( -2pi*i*f*k ) )
                    // Where f = sampling frequency ( = bin # / N )
                    // and with our FFT and DFT t_k = k (sample frequency is the index of the array)
                    // 
                    // At this stage in the FFT algorithm, N is actually the length of each
                    // independent butterfly chain (seperation).
                    //
                    // For solving euler's identity, we use real and imaginary components such that:
                    // e^( aik ) = cos(aik) + isin(aik), where a is a constant coefficient, i is sqrt(-1)
                    // and k is the parameter.
                    // 
                    // So then e^( -2pi*i*f*k ) = cos( -2pi*f*k ) - sin( -2pi*f*k )
                    // Note that sin is negated since we travel counter clockwise (-) around the complex plane
                    //
                    // There are some possible optimisations here, but for sake of understanding
                    // algorithm and math we will avoid them
                    //
                    // Finally, the DFT of a value can be calculated as follows:
                    // DFT(k) = g_r[t_k] * e^( -2pi*i*f*k ) - g_i[t_k] * e^( 2pi*i*f*k )
                    // Note the positive 2pi in the second euler identity.
                    //
                    // This expands to the following using Euler's formulas:
                    // DFT(k) = g_r[t_k] * [ cos( 2pi*f*k ) + sin ( 2pi*f*k ) ]
                    // - g_i[t_k] * [ -cos( 2pi*f*k ) + sin( 2pi*f*k ) ]
                    //
                    // Expanded:
                    // = [R] g_r[t_k] * cos( 2pi*f*k ) + [I] g_r[t_k] * sin( 2pi*f*k )
                    // + [I] g_i[t_k] * cos( 2pi*f*k ) - [R] g_i[t_k] * sin( 2pi*f*k )
                    //
                    // With the following real and imaginary components:
                    // R = g_r[t_k] * cos( 2pi*f*k ) - g_i[t_k] * sin( 2pi*f*k )
                    // I = g_r[t_k] * sin( 2pi*f*k ) + g_i[t_k] * cos( 2pi*f*k )
                    //
                    // So we end up with the following real and imaginary values for lo, 
                    // where cos and sin are precalculated twiddle factors:
                    let valLoR = inReal[lo] * twiddleReal - inImaginary[lo] * twiddleImaginary;
                    let valLoI = inReal[lo] * twiddleImaginary + inImaginary[lo] * twiddleReal;

                    // 2. Set lo = hi - DFT(lo)
                    inReal[lo] = inReal[hi] - valLoR;
                    inImaginary[lo] = inImaginary[hi] - valLoI;
                    // 3. Set hi = hi + DFT(lo)
                    inReal[hi] += valLoR;
                    inImaginary[hi] += valLoI;
                }

                // Rotate twiddle 1/seperation way, counter-clockwise (-), around euler's circle on complex plane
                // Recall clockwise: 
                //      RealValue =  RealValue*cos(2pif) + ImaginaryValue*sin(2pif)
                // ImaginaryValue = -RealValue*sin(2pif) + ImaginaryValue*cos(2pif)
                // And counter-clockwise:
                //      RealValue = RealValue*cos(2pif) - ImaginaryValue*sin(2pif)
                // ImaginaryValue = RealValue*sin(2pif) + ImaginaryValue*cos(2pif)
                let t = twiddleReal;
                twiddleReal = t * Math.cos(TWOPI/separation) - direction * twiddleImaginary * Math.sin(TWOPI/separation);
                twiddleImaginary = direction * t * Math.sin(TWOPI/separation) + twiddleImaginary * Math.cos(TWOPI/separation);
            }
        }
    }
    else if (mode == 1){
        // Radix 2
        // Log2 size stages (see diagram for how many butterfly stages there are)
        for (let stage = 1; stage <= logSize; stage++) {
            // Separation between butterfly chains = total length of each butterfly chain
            let separation = Math.pow(2, stage);
            // The width of each butterfly chain = distance between lo and hi in DFT
            let width = separation / 2;
            // For each offset of butterfly chain
            for (let i = 0; i < width; i++) {
                // ... in each butterfly chain
                for (let currentButterfly = i; currentButterfly < size; currentButterfly += separation) {

                    // Get hi and lo ptrs for DFT
                    let hi = currentButterfly;
                    let lo = hi + width;

                    // Algorithm states to take DFT of lo and add it to hi, and negate
                    // the value at lo and add the hi

                    // 1. Get value of lo (rel and imaginary)
                    // Recall DFT = sum_k=0^N ( g(t_k) * e^( -2pi*i*f*k ) )
                    // Where f = sampling frequency ( = bin # / N )
                    // and with our FFT and DFT t_k = k (sample frequency is the index of the array)
                    // 
                    // At this stage in the FFT algorithm, N is actually the length of each
                    // independent butterfly chain (seperation).
                    //
                    // For solving euler's identity, we use real and imaginary components such that:
                    // e^( aik ) = cos(aik) + isin(aik), where a is a constant coefficient, i is sqrt(-1)
                    // and k is the parameter.
                    // 
                    // So then e^( -2pi*i*f*k ) = cos( -2pi*f*k ) - sin( -2pi*f*k )
                    // Note that sin is negated since we travel counter clockwise (-) around the complex plane
                    //
                    // There are some possible optimisations here, but for sake of understanding
                    // algorithm and math we will avoid them
                    //
                    // Finally, the DFT of a value can be calculated as follows:
                    // DFT(k) = g_r[t_k] * e^( -2pi*i*f*k ) - g_i[t_k] * e^( 2pi*i*f*k )
                    // Note the positive 2pi in the second euler identity.
                    // 
                    // from https://www.analog.com/media/en/technical-documentation/dsp-book/dsp_book_Ch31.pdf
                    //
                    // This expands to the following using Euler's formulas:
                    // DFT(k) = g_r[t_k] * [ cos( 2pi*f*k ) + sin ( 2pi*f*k ) ]
                    // - g_i[t_k] * [ -cos( 2pi*f*k ) + sin( 2pi*f*k ) ]
                    //
                    // Expanded:
                    // = [R] g_r[t_k] * cos( 2pi*f*k ) + [I] g_r[t_k] * sin( 2pi*f*k )
                    // + [I] g_i[t_k] * cos( 2pi*f*k ) - [R] g_i[t_k] * sin( 2pi*f*k )
                    //
                    // With the following real and imaginary components:
                    // R = g_r[t_k] * cos( 2pi*f*k ) - g_i[t_k] * sin( 2pi*f*k )
                    // I = g_r[t_k] * sin( 2pi*f*k ) + g_i[t_k] * cos( 2pi*f*k )
                    //
                    // Also note the frequency is = i (the offset of butterfly chain) / separation
                    // because for each stage the sample size is effectively = separation 
                    // ( that is, separation = the length from one butterfly to the next = # of elements in
                    // each one )
                    let f = i / separation;
                    //
                    // So we end up with the following real and imaginary values for lo:
                    let valLoR = inReal[lo] * Math.cos(TWOPI*f) - inImaginary[lo] * direction * Math.sin(TWOPI*f);
                    let valLoI = inReal[lo] * direction * Math.sin(TWOPI*f) + inImaginary[lo] * Math.cos(TWOPI*f);

                    // 2. Set lo = hi - DFT(lo)
                    inReal[lo] = inReal[hi] - valLoR;
                    inImaginary[lo] = inImaginary[hi] - valLoI;
                    // 3. Set hi = hi + DFT(lo)
                    inReal[hi] += valLoR;
                    inImaginary[hi] += valLoI;
                }
            }
        }
    }
    else if (mode == 2) {

        // Radix 2
        // Log2 size stages (see diagram for how many butterfly stages there are)
        for (let stage = 1; stage <= logSize; stage++) {

            // distance between the start of each butterfly chain
            let butterflySeperation = Math.pow(2,stage);

            // distance between lo and hi elements in butterfly
            let butterflyWidth = butterflySeperation / 2; 

            // For each butterfly separation
            for (let j = 0; j < size; j += butterflySeperation) {

                // Perform synthesis on offset (i-j) and offset + width (hi and lo) for each butterfly chain
                for (let i = j; i < j + butterflyWidth; i++) {

                    // get frequency, set coefficient
                    // Note that sampling size is the butterfly seperation
                    let f = (i-j) / butterflySeperation;
                    let a = TWOPI * f;

                    // Indices for hi and low values (low is actually further up)
                    let idxHi = i;
                    let idxLo = idxHi + butterflyWidth;

                    // Get all time values, real and imaginary
                    let gtkrHi = inReal[idxHi];
                    let gtkiHi = inImaginary[idxHi];
                    let gtkrLo = inReal[idxLo];
                    let gtkiLo = inImaginary[idxLo];
                    
                    // Recall the below is the exact same as adding all vals from DFT
                    let ValRLo = gtkrLo * Math.cos(a) - direction * gtkiLo * Math.sin(a);
                    let ValILo = direction * gtkrLo * Math.sin(a) + gtkiLo * Math.cos(a);
                    // Nature of the Radix 2 algorithm means the hi val is equal to itself (see diagram)
                    let ValRHi = gtkrHi;
                    let ValIHi = gtkiHi;

                    // Now add the lo to the high and TODO
                    inReal[idxHi] = ValRLo + ValRHi;
                    inImaginary[idxHi] = ValILo + ValIHi;
                    inReal[idxLo] = ValRHi - ValRLo;
                    inImaginary[idxLo] = ValIHi - ValILo;
                }
            }
        }
    }
    else if (mode == 3) {
        let tr, ti;
        let len = size;
        let arr = [];
        let M_PI = Math.PI;
        let inv = direction == -1;
        for (let i = 0; i < size; i++) {
            arr[i] = inReal[i];
            arr[i+len] = inImaginary[i];
        }
        for (let i = 1; i <= Math.log2(len); i++) {
            let le = Math.pow(2, i);
            let le2 = le / 2;
            let ur = 1;
            let ui = 0;
            let sr = Math.cos(-M_PI / le2);
            let si = direction * -Math.sin(-M_PI / le2);
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
    }
    // It was in place so set references
    outReal = inReal;
    outImaginary = inImaginary;


    console.log("FFT Mode " + mode + " completed in " + (window.performance.now() - time) + " milliseconds!")

    // If inverse we have to divide by the size of input (our reference)
    if (direction == -1) {
        for (let i = 0; i < size; i++) {
            outReal[i] /= size;
            outImaginary[i] /= size;
        }
    }

    // Now we apply to out table
    let result = {};
    result.real = outReal;
    result.imaginary = outImaginary;
    return result;
}