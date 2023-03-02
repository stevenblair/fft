package fft

import "math"

const (
	w61 = 0.5 - 0.86602540378443864676372317075293618347140262690519031402i  // exp(-2πi*1/6)
	w62 = -0.5 - 0.86602540378443864676372317075293618347140262690519031402i // exp(-2πi*2/6)
	w63 = -1                                                                 // exp(-2πi*3/6)
	w64 = -0.5 + 0.86602540378443864676372317075293618347140262690519031402i // exp(-2πi*4/6)
	w65 = 0.5 + 0.86602540378443864676372317075293618347140262690519031402i  // exp(-2πi*5/6)
)

func (f *FFTData) radix6(x []complex128, s int) []complex128 {
	const r = 6
	n := len(x)

	// Copy input into new storage.
	copy(f.y, x)

	// Reorder input using base-r digit reversal permutation.
	for i, j := 0, 0; i < n-1; i++ {
		if i < j {
			f.y[i], f.y[j] = f.y[j], f.y[i]
		}
		k := (r - 1) * n / r
		for k <= j {
			j -= k
			k /= r
		}
		j += k / (r - 1)
	}

	for m := r; m <= n; m *= r {
		// Calculate twiddle factor.
		wim, wre := math.Sincos(-2 * math.Pi * float64(s) / float64(m))
		w := complex(wre, wim)

		mr := m / r
		for i, wi := 0, complex(1, 0); i < mr; i++ {
			for j := 0; j < n; j += m {
				// Retrieve subset of points.
				t0, t1, t2, t3, t4, t5 := f.y[i+j], f.y[i+j+mr], f.y[i+j+2*mr], f.y[i+j+3*mr], f.y[i+j+4*mr], f.y[i+j+5*mr]

				// Apply twiddle factors w**(i+k) for 1 ≤ k < r.
				t1 *= wi
				t2 *= wi * wi
				t3 *= wi * wi * wi
				t4 *= wi * wi * wi * wi
				t5 *= wi * wi * wi * wi * wi

				// Transform points using r-point DFT.
				f.y[i+j] += t1 + t2 + t3 + t4 + t5
				f.y[i+j+3*mr] = t0 - t1 + t2 - t3 + t4 - t5
				if s > 0 {
					f.y[i+j+mr], f.y[i+j+2*mr], f.y[i+j+4*mr], f.y[i+j+5*mr] =
						t0+t1*w61+t2*w62-t3+t4*w64+t5*w65,
						t0+t1*w62+t2*w64+t3+t4*w62+t5*w64,
						t0+t1*w64+t2*w62+t3+t4*w64+t5*w62,
						t0+t1*w65+t2*w64-t3+t4*w62+t5*w61
				} else {
					// 1/w61 = w65, etc.
					f.y[i+j+mr], f.y[i+j+2*mr], f.y[i+j+4*mr], f.y[i+j+5*mr] =
						t0+t1*w65+t2*w64-t3+t4*w62+t5*w61,
						t0+t1*w64+t2*w62+t3+t4*w64+t5*w62,
						t0+t1*w62+t2*w64+t3+t4*w62+t5*w64,
						t0+t1*w61+t2*w62-t3+t4*w64+t5*w65
				}
			}
			wi *= w
		}
	}
	return f.y
}
