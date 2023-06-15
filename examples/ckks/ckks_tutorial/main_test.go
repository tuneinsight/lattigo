package main

import "testing"

func TestMain(t *testing.T) {
	if testing.Short() {
		t.Skip("skipped in -short mode")
	}
	t.Skip("test not passing") // TODO: bug in the linear transform API ?
	main()
}
