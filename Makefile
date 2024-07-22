.DEFAULT_GOAL := test

.PHONY: test_gotest
test_gotest:
	go clean -testcache
	go test -timeout=0 ./...

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
	
	@STATICCHECKOUT=$$(staticcheck -go 1.22 -checks all ./...); \
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
test: test_gotest

.PHONY: ci_test
ci_test: static_check test_gotest

EXECUTABLES = goimports staticcheck
.PHONY: get_tools
get_tools:
	go install golang.org/x/tools/cmd/goimports@latest
	go install honnef.co/go/tools/cmd/staticcheck@2023.1.7

.PHONY: check_tools
check_tools:
	@$(foreach exec,$(EXECUTABLES),\
		$(if $(shell which $(exec)),true,$(error "$(exec) not found in PATH, consider running `make get_tools`.")))