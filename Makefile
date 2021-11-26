.DEFAULT_GOAL := test

Coding/bin/Makefile.base:
	git clone https://github.com/dedis/Coding

.PHONY: test_examples
test_examples:
	@echo Running the examples
	go run ./examples/ring/vOLE -short > /dev/null
	go run ./examples/bfv > /dev/null
	go run ./examples/ckks/euler > /dev/null
	go run ./examples/ckks/sigmoid > /dev/null
	go run ./examples/rlwe/lwe_bridge > /dev/null
	go run ./examples/dbfv/pir &> /dev/null
	go run ./examples/dbfv/psi &> /dev/null
	@echo ok
	@echo Building resources-heavy examples
	go build -o /dev/null ./examples/ckks/bootstrapping
	go build -o /dev/null ./examples/ckks/advanced/rlwe_lwe_bridge_LHHMQ20
	@echo ok

.PHONY: test_gotest
test_gotest:
	go test -v -timeout=0 ./utils ./ring ./bfv ./ckks ./dbfv ./dckks
	go test -v -timeout=0 ./ckks/advanced
	go test -v -timeout=0 ./ckks/bootstrapping -test-bootstrapping -short

.PHONY: test
test: test_fmt test_gotest test_examples

.PHONY: ci_test
ci_test: test_fmt test_lint test_gotest test_examples

%: force Coding/bin/Makefile.base
	@$(MAKE) -f Coding/bin/Makefile.base $@
force: ;