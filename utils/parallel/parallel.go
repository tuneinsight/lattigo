// Package parallel provides a process-wide bound on the number of goroutines
// used to evaluate homomorphic circuits (linear transformations, polynomials,
// bootstrapping, RNS-limb-wise ring arithmetic) and a helper to distribute
// independent work.
//
// Evaluators are safe for concurrent use since v6.2.0 (all temporary buffers
// are drawn from sync.Pool-backed pools), which makes this distribution
// possible. The bound defaults to runtime.GOMAXPROCS(0); applications that
// already parallelize at a coarser level can call SetParallelism(1) to restore
// fully sequential evaluation.
//
// Parallel regions may be nested (e.g. the giant steps of a linear
// transformation each run RNS-limb-parallel NTTs). To avoid oversubscription,
// helper goroutines are drawn from a global token pool of size
// Parallelism()-1: the calling goroutine always works itself and only borrows
// helpers that are idle, so nested regions degrade gracefully to sequential
// execution when all cores are busy.
package parallel

import (
	"runtime"
	"sync"
	"sync/atomic"
)

type tokenPool struct {
	size   int
	tokens chan struct{}
}

func newTokenPool(n int) *tokenPool {
	if n < 1 {
		n = 1
	}
	p := &tokenPool{size: n, tokens: make(chan struct{}, n-1)}
	for i := 0; i < n-1; i++ {
		p.tokens <- struct{}{}
	}
	return p
}

var pool atomic.Pointer[tokenPool]

func init() {
	pool.Store(newTokenPool(runtime.GOMAXPROCS(0)))
}

// SetParallelism sets the maximum number of goroutines used concurrently by
// parallel regions. Values < 1 are clamped to 1 (sequential).
func SetParallelism(n int) {
	pool.Store(newTokenPool(n))
}

// Parallelism returns the maximum number of goroutines used concurrently by
// parallel regions.
func Parallelism() int {
	return pool.Load().size
}

// For runs f(0..n-1), distributing the indices over the calling goroutine
// and as many idle helper goroutines as available (at most n-1), and returns
// the first error encountered. f must be safe to call concurrently for
// distinct indices.
//
// Helpers are borrowed when the region starts and, while unclaimed indices
// remain, whenever a worker claims an index: tokens returned by regions that
// finished in the meantime join regions still running. Otherwise a region
// started while the pool was busy (e.g. the largest of several concurrent
// linear transformations) would keep its few helpers until it finished.
func For(n int, f func(i int) error) error {
	p := pool.Load()

	var next atomic.Int64
	var workers atomic.Int64 // calling goroutine and helpers
	var failed atomic.Bool
	var wg sync.WaitGroup

	var mu sync.Mutex
	var firstErr error
	var firstPanic any

	// tryBorrow takes an idle token without blocking if fewer workers than
	// unclaimed indices are running.
	tryBorrow := func() bool {
		if workers.Load() >= int64(n)-next.Load() {
			return false
		}
		select {
		case <-p.tokens:
			workers.Add(1)
			return true
		default:
			return false
		}
	}

	var spawn func()
	work := func() {
		// A panic in a helper goroutine could not be recovered by the
		// caller and would terminate the process. Every worker (including
		// the calling goroutine) therefore records its panic, and the first
		// one is re-raised in the caller once all workers have stopped.
		defer func() {
			if r := recover(); r != nil {
				mu.Lock()
				if firstPanic == nil {
					firstPanic = r
				}
				mu.Unlock()
				failed.Store(true)
			}
		}()
		for !failed.Load() {
			i := int(next.Add(1) - 1)
			if i >= n {
				return
			}
			if tryBorrow() {
				spawn()
			}
			if err := f(i); err != nil {
				mu.Lock()
				if firstErr == nil {
					firstErr = err
				}
				mu.Unlock()
				failed.Store(true)
				return
			}
		}
	}
	spawn = func() {
		wg.Add(1)
		go func() {
			defer func() {
				p.tokens <- struct{}{}
				wg.Done()
			}()
			work()
		}()
	}

	// Borrow idle helpers without blocking.
	workers.Store(1)
	for workers.Load() < int64(n) && tryBorrow() {
		spawn()
	}

	work()
	wg.Wait()

	if firstPanic != nil {
		panic(firstPanic)
	}
	return firstErr
}

// ForRange splits [0, n) into contiguous chunks of at least minChunk
// elements and runs f(start, end) on each chunk in parallel.
func ForRange(n, minChunk int, f func(start, end int)) {
	chunks := min(Parallelism(), max(n/max(minChunk, 1), 1))
	if chunks <= 1 {
		f(0, n)
		return
	}
	size := (n + chunks - 1) / chunks
	_ = For(chunks, func(c int) error {
		start := c * size
		f(start, min(start+size, n))
		return nil
	})
}
