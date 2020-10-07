.DEFAULT_GOAL := local

Coding/bin/Makefile.base:
	git clone https://github.com/dedis/Coding
include Coding/bin/Makefile.base

.PHONY: test_examples
test_examples:
	@echo Running the examples
	go run ./examples/bfv > /dev/null
	go run ./examples/ckks/euler > /dev/null
	go run ./examples/ckks/sigmoid > /dev/null
	go run ./examples/dbfv/pir &> /dev/null
	go run ./examples/dbfv/psi &> /dev/null
	@echo ok
	@echo Building resources-heavy examples
	go build -o /dev/null ./examples/ckks/bootstrapping
	@echo ok

.PHONY: test_gotest
	go test -v -short -p=1 ./... -timeout=0

.PHONY: test
test: test_fmt test_gotest test_examples

.PHONY: local
local: test_fmt test_lint test_gotest test_examples
