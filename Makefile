.PHONY: install test

install:
	pip install -e .[dev]

test:
	PYTHONPATH=src python -m unittest discover -s tests
