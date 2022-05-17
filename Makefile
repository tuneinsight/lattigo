.DEFAULT_GOAL := test

.PHONY: test_gotest
test_gotest:
	go test -v -timeout=0 ./utils ./ring ./bfv ./ckks ./dbfv ./dckks
	go test -v -timeout=0 ./ckks/advanced
	go test -v -timeout=0 ./ckks/bootstrapping -test-bootstrapping -short

.PHONY: test_examples
test_examples:
	@echo Running the examples
	go run ./examples/ring/vOLE -short > /dev/null
	go run ./examples/bfv > /dev/null
	go run ./examples/ckks/euler > /dev/null
	go run ./examples/ckks/polyeval > /dev/null
	go run ./examples/rlwe/lwe_bridge > /dev/null
	go run ./examples/dbfv/pir &> /dev/null
	go run ./examples/dbfv/psi &> /dev/null
	@echo ok
	@echo Building resources-heavy examples
	go build -o /dev/null ./examples/ckks/bootstrapping
	go build -o /dev/null ./examples/ckks/advanced/rlwe_lwe_bridge_LHHMQ20
	@echo ok

.PHONY: static_checks
static_checks:
	@echo Checking correct formatting of files
	out=`go fmt ./...`; echo "$$out"; [ -z "$$out" ]
	go vet ./...
	out=`golint ./...`; echo "$$out"; [ -z "$$out" ]
	staticcheck -go 1.17 ./...
	go mod tidy
	out=`git status --porcelain`; echo "$$out"; [ -z "$$out" ]

.PHONY: test
test: test_gotest test_examples

.PHONY: ci_test
ci_test: static_checks test_gotest test_examples

.PHONY: get_tools
get_tools:
	go install golang.org/x/lint/golint@latest
	go install honnef.co/go/tools/cmd/staticcheck@2022.1.1