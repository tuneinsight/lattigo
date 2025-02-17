.DEFAULT_GOAL := test

.PHONY: test_gotest
test_gotest:
	go clean -testcache
	go test -timeout=0 ./...

.PHONY: checks
checks: check_tools
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
	
	@STATICCHECKOUT=$$(staticcheck -go 1.23 -checks all ./...); \
	if [ -z "$$STATICCHECKOUT" ]; then\
        echo "staticcheck: OK";\
	else \
		echo "staticcheck: problems in files:";\
		echo "$$STATICCHECKOUT";\
		false;\
    fi

	@GOVULNCHECKOUT=$$(govulncheck ./...); \
	if echo "$$GOVULNCHECKOUT" | grep -q "No vulnerabilities found"; then\
		echo "govulncheck: OK";\
    else \
		echo "govulncheck:" >&2;\
		echo "$$GOVULNCHECKOUT" >&2;\
		false;\
	fi

	@GOSECOUT=$$(gosec -quiet ./...); \
	if [ -z "$$GOSECOUT" ]; then\
		echo "gosec: OK";\
	else \
		echo "gosec: problems in files:";\
		echo "$$GOSECOUT";\
		false;\
	fi
	
	@echo Checking all local changes are committed
	go mod tidy
	out=`git status --porcelain`; echo "$$out"; [ -z "$$out" ]

.PHONY: test
test: test_gotest

.PHONY: ci_test
ci_test: checks test_gotest

EXECUTABLES = goimports staticcheck govulncheck gosec
.PHONY: get_tools
get_tools:
	go install golang.org/x/tools/cmd/goimports@latest
	go install honnef.co/go/tools/cmd/staticcheck@2024.1.1
	go install golang.org/x/vuln/cmd/govulncheck@latest
	go install github.com/securego/gosec/v2/cmd/gosec@latest

.PHONY: check_tools
check_tools:
	@$(foreach exec,$(EXECUTABLES),\
		$(if $(shell which $(exec)),true,$(error "$(exec) not found in PATH, consider running `make get_tools`.")))