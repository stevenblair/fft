package fft

import "math"

// Bluestein algorithm for prime-n sized FFTs, i.e. pad input to a size m = 2**p ≥ 2*n-1
// for some integer p and use the circular convolution theorem to turn the size n DFT
// into two size m DFTs and a size m inverse DFT.

func (f *FFTData) bluestein(x []complex128) []complex128 {
	copy(f.y2, x)

	a0 := math.Pi / float64(f.n)
	f.w[0] = 1
	for i := 1; i < f.n; i++ {
		s, c := math.Sincos(a0 * float64(i*i))
		f.w[i] = complex(c, s)
		f.w[f.m-i] = complex(c, s)
		f.y2[i] *= complex(c, -s)
	}

	f.s1.y = f.s1.stockham(f.y2, 1)
	for i, ww := range f.s2.stockham(f.w, 1) {
		f.s1.y[i] *= ww
	}
	f.s1.y = f.s1.stockham(f.s1.y, -1)

	for i := 0; i < f.n; i++ {
		f.s1.y[i] *= complex(real(f.w[i])/float64(f.m), -imag(f.w[i])/float64(f.m))
	}

	return f.s1.y[:f.n]
}

func (f *FFTData) bluesteini(x []complex128) []complex128 {
	for i, xi := range x {
		f.y[i] = complex(real(xi), -imag(xi))
	}
	f.y = f.bluestein(f.y)
	for i, yi := range f.y {
		f.y[i] = complex(real(yi)/float64(f.n), -imag(yi)/float64(f.n))
	}
	return f.y
}
