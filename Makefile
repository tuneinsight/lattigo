.DEFAULT_GOAL := local

Coding/bin/Makefile.base:
	git clone https://github.com/dedis/Coding
include Coding/bin/Makefile.base

.PHONY: test_examples
test_examples:
	go run ./examples/bfv/examples_bfv.go
	go run ./examples/ckks/examples_ckks.go

.PHONY: test_local
test_local:
	go test -v ./...

.PHONY: local
local: test_fmt test_lint test_local test_examples
