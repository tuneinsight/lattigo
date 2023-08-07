package main

import "testing"

func TestMain(t *testing.T) {
	if testing.Short() {
		t.Skip("skipped in -short mode")
	}
	main()
}
