package fft

import "math"

// Bluestein algorithm for prime-n sized FFTs, i.e. pad input to a size m = 2**p ≥ 2*n-1
// for some integer p and use the circular convolution theorem to turn the size n DFT
// into two size m DFTs and a size m inverse DFT.

func (f *FFTData) bluestein(x []complex128) []complex128 {
	n := len(x)

	m := 1 << uint(math.Ilogb(float64(2*n-1)))
	if m < 2*n-1 {
		m <<= 1
	}

	// w := make([]complex128, m)
	// y := make([]complex128, m)
	copy(f.y, x)

	a0 := math.Pi / float64(n)
	f.w[0] = 1
	for i := 1; i < n; i++ {
		s, c := math.Sincos(a0 * float64(i*i))
		f.w[i] = complex(c, s)
		f.w[m-i] = complex(c, s)
		f.y[i] *= complex(c, -s)
	}

	f.y1.y = f.y1.stockham(f.y, 1)
	for i, ww := range f.y2.stockham(f.w, 1) {
		f.y2.y[i] *= ww
	}
	f.y3.y = f.y3.stockham(f.y2.y, -1)

	for i := 0; i < n; i++ {
		f.y3.y[i] *= complex(real(f.w[i])/float64(m), -imag(f.w[i])/float64(m))
	}

	return f.y3.y[:n]
}

func (f *FFTData) bluesteini(x []complex128) []complex128 {
	n := len(x)
	// y := make([]complex128, n)
	for i, xi := range x {
		f.y[i] = complex(real(xi), -imag(xi))
	}
	f.y = f.bluestein(f.y)
	for i, yi := range f.y {
		f.y[i] = complex(real(yi)/float64(n), -imag(yi)/float64(n))
	}
	return f.y
}
