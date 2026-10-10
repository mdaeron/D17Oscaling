main:
	clear; uv run python src/D17Oscaling/__init__.py

build-example-data:
	uv run examples/build-example-data.py
