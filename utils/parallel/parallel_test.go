package parallel

import (
	"errors"
	"sync/atomic"
	"testing"
	"time"
)

func TestForRunsEveryIndexOnce(t *testing.T) {
	defer SetParallelism(Parallelism())
	SetParallelism(8)
	var hits [1000]atomic.Int32
	if err := For(len(hits), func(i int) error { hits[i].Add(1); return nil }); err != nil {
		t.Fatal(err)
	}
	for i := range hits {
		if hits[i].Load() != 1 {
			t.Fatalf("index %d ran %d times", i, hits[i].Load())
		}
	}
}

func TestForReturnsError(t *testing.T) {
	defer SetParallelism(Parallelism())
	SetParallelism(8)
	want := errors.New("boom")
	if err := For(100, func(i int) error {
		if i == 57 {
			return want
		}
		return nil
	}); !errors.Is(err, want) {
		t.Fatalf("got %v, want %v", err, want)
	}
}

func TestForReraisesHelperPanicInCaller(t *testing.T) {
	defer SetParallelism(Parallelism())
	SetParallelism(8)
	defer func() {
		if r := recover(); r != "boom" {
			t.Fatalf("recovered %v, want boom", r)
		}
		// All tokens must have been returned.
		if got := len(pool.Load().tokens); got != 7 {
			t.Fatalf("%d of 7 tokens returned", got)
		}
	}()
	_ = For(100, func(i int) error {
		if i == 99 { // most likely executed by a helper goroutine
			panic("boom")
		}
		return nil
	})
	t.Fatal("For did not panic")
}

func TestNestedForDoesNotDeadlock(t *testing.T) {
	defer SetParallelism(Parallelism())
	SetParallelism(4)
	var sum atomic.Int64
	_ = For(16, func(i int) error {
		return For(16, func(j int) error { sum.Add(1); return nil })
	})
	if sum.Load() != 256 {
		t.Fatalf("sum = %d", sum.Load())
	}
}

func TestForBorrowsReleasedTokens(t *testing.T) {
	defer SetParallelism(Parallelism())
	SetParallelism(4)

	// Hold all tokens, as concurrent regions would.
	p := pool.Load()
	for i := 0; i < 3; i++ {
		<-p.tokens
	}

	var active, peak atomic.Int32
	released := make(chan struct{})
	err := For(8, func(i int) error {
		if a := active.Add(1); a >= 2 {
			peak.Store(a)
		}
		defer active.Add(-1)
		if i == 0 {
			// The region started without helpers; release the tokens.
			for j := 0; j < 3; j++ {
				p.tokens <- struct{}{}
			}
			close(released)
			return nil
		}
		<-released
		// A helper borrowed after the release must join this region.
		deadline := time.Now().Add(5 * time.Second)
		for peak.Load() < 2 {
			if time.Now().After(deadline) {
				return errors.New("no helper joined the running region")
			}
			time.Sleep(time.Millisecond)
		}
		return nil
	})
	if err != nil {
		t.Fatal(err)
	}
	if got := len(p.tokens); got != 3 {
		t.Fatalf("%d of 3 tokens returned", got)
	}
}
