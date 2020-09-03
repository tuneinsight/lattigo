EXCLUDE_LINT = "_test.go"

test_fmt:
	@echo Checking correct formatting of files
	@{ \
		files=$$( go fmt ./... ); \
		if [ -n "$$files" ]; then \
		echo "Files not properly formatted: $$files"; \
		exit 1; \
		fi; \
		if ! go vet ./...; then \
		exit 1; \
		fi \
	}

test_lint:
	@echo Checking linting of files
	@{ \
		go install golang.org/x/lint/golint; \
		el=$(EXCLUDE_LINT); \
		lintfiles=$$( golint ./... | egrep -v "$$el" ); \
		if [ -n "$$lintfiles" ]; then \
		echo "Lint errors:"; \
		echo "$$lintfiles"; \
		exit 1; \
		fi \
	}

test_local:
	go test -v -short -p=1 ./... -timeout=0
	@echo Running the exemples
	go run ./examples/bfv/examples_bfv.go > /dev/null
	go run ./examples/ckks/euler/euler.go > /dev/null
	go run ./examples/ckks/sigmoid/sigmoid.go > /dev/null
	go run ./examples/dbfv/pir/pir.go &> /dev/null
	go run ./examples/dbfv/psi/psi.go &> /dev/null
	@echo ok

test: test_fmt test_local

local: test_fmt test_lint test_local
