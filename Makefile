.DEFAULT_GOAL := test

Coding/bin/Makefile.base:
	git clone https://github.com/dedis/Coding

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
test_gotest:
	go test -v -timeout=0 ./utils ./ring ./bfv ./dbfv ./dckks
	go test -v -timeout=0 ./ckks -test-bootstrapping

.PHONY: test
test: test_fmt test_gotest test_examples

.PHONY: ci_test
ci_test: test_fmt test_lint test_gotest test_examples

%: force Coding/bin/Makefile.base
	@$(MAKE) -f Coding/bin/Makefile.base $@
force: ;