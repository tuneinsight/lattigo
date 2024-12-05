help: ## Display this help screen
	@grep -h \
		-E '^[a-zA-Z_0-9-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

test_gotest: ## test_gotest
	go clean -testcache
	go test -timeout=0 ./...

checks: check_tools ## checks
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

	@GOVULNCHECKOUT=$$(govulncheck ./...); \
	if echo "$$GOVULNCHECKOUT" | grep -q "No vulnerabilities found"; then\
		echo "govulncheck: OK";\
    else \
		echo "govulncheck:" >&2;\
		echo "$$GOVULNCHECKOUT" >&2;\
		false;\
	fi

# gosec rule G115: Is exluded because there are int->uin64 conversions
# and the rule currently contains false positives
	@GOSECOUT=$$(gosec -quiet -exclude=G115 ./...); \
	if [ -z "$$GOSECOUT" ]; then\
		echo "gosec: OK (excluding G115)";\
	else \
		echo "gosec: problems in files:";\
		echo "$$GOSECOUT";\
		false;\
	fi
	
	@echo Checking all local changes are committed
	go mod tidy
	out=`git status --porcelain`; echo "$$out"; [ -z "$$out" ]

test: test_gotest ## test

ci_test: checks test_gotest ## ci_test

EXECUTABLES = goimports staticcheck govulncheck gosec
get_tools: ## get tools
	go install golang.org/x/tools/cmd/goimports@latest
	go install honnef.co/go/tools/cmd/staticcheck@2023.1.7
	go install golang.org/x/vuln/cmd/govulncheck@latest
	go install github.com/securego/gosec/v2/cmd/gosec@latest


check_tools: ## check tools
	@$(foreach exec,$(EXECUTABLES),\
		$(if $(shell which $(exec)),true,$(error "$(exec) not found in PATH, consider running `make get_tools`.")))


.PHONY: test_gotest checks  check_tools get_tools ci_test test help

