# Adds file annotations to Github Actions (only useful on CI)
GITHUB_ACTIONS_FORMATTING=0
ifeq ($(GITHUB_ACTIONS_FORMATTING), 1)
	FLAKE8_FORMAT=--format='::error file=%(path)s,line=%(row)d,col=%(col)d,title=%(code)s::%(path)s:%(row)d:%(col)d: %(code)s %(text)s'
else
	FLAKE8_FORMAT=
endif

doc:
	@(cd docs/ && make html)

lint:
	@echo "    Linting FUSE codebase"
	@python -m flake8 $(FLAKE8_FORMAT) fuse
	@echo "    Linting FUSE test suite"
	@python -m flake8 $(FLAKE8_FORMAT) test

test_examples:
	@echo "    Running examples"
	@python -m pytest test/test_2d_examples_docs.py
	@python -m pytest test/test_3d_examples_docs.py

tests:
	@echo "    Running all tests"
	@python -m coverage run -p -m pytest -rx test

coverage:
	@python -m coverage combine
	@python -m coverage report -m --skip-covered
	@python -m coverage json

test_cells:
	@echo "    Running all cell comparison tests"
	@firedrake-clean
	@python -m pytest -rPx --run-cleared test/test_cells.py::test_ref_els[expect0]
	@firedrake-clean
	@python -m pytest -rPx --run-cleared test/test_cells.py::test_ref_els[expect1]

prepush: lint tests doc