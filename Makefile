.DEFAULT_GOAL := test

.PHONY: test_gotest
test_gotest:
	go test -timeout=0 ./...

.PHONY: test_examples
test_examples:
	@echo Running the examples
	go run ./examples/ring/vOLE -short > /dev/null
	go run ./examples/rgsw > /dev/null
	go run ./examples/bfv > /dev/null
	go run ./examples/ckks/bootstrapping -short > /dev/null
	go run ./examples/ckks/advanced/lut -short > /dev/null
	go run ./examples/ckks/euler > /dev/null
	go run ./examples/ckks/polyeval > /dev/null
	go run ./examples/dbfv/pir &> /dev/null
	go run ./examples/dbfv/psi &> /dev/null
	@echo ok

.PHONY: static_check
static_check: check_tools
	@echo Checking correct formatting of files
	
	@FMTOUT=$$(go fmt ./...); \
	if [ -z $$FMTOUT ]; then\
        echo "go fmt: OK";\
	else \
		echo "go fmt: problems in files:";\
		echo $$FMTOUT;\
		false;\
    fi

	@if GOVETOUT=$$(go vet ./... 2>&1); then\
        echo "go vet: OK";\
	else \
		echo "go vet: problems in files:";\
		echo "$$GOVETOUT";\
		false;\
    fi

	@GOIMPORTSOUT=$$(goimports -l .); \
	if [ -z "$$GOIMPORTSOUT" ]; then\
        echo "goimports: OK";\
	else \
		echo "goimports: problems in files:";\
		echo "$$GOIMPORTSOUT";\
		false;\
    fi
	
	@STATICCHECKOUT=$$(staticcheck -go 1.20 -checks all ./...); \
	if [ -z "$$STATICCHECKOUT" ]; then\
        echo "staticcheck: OK";\
	else \
		echo "staticcheck: problems in files:";\
		echo "$$STATICCHECKOUT";\
		false;\
    fi
	
	@echo Checking all local changes are committed
	go mod tidy
	out=`git status --porcelain`; echo "$$out"; [ -z "$$out" ]

.PHONY: test
test: test_gotest test_examples

.PHONY: ci_test
ci_test: static_check test_gotest test_examples

EXECUTABLES = goimports staticcheck
.PHONY: get_tools
get_tools:
	go install golang.org/x/tools/cmd/goimports@latest
	go install honnef.co/go/tools/cmd/staticcheck@latest

.PHONY: check_tools
check_tools:
	@$(foreach exec,$(EXECUTABLES),\
		$(if $(shell which $(exec)),true,$(error "$(exec) not found in PATH, consider running `make get_tools`.")))