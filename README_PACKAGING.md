Packaging and local install

Build a wheel and install locally for testing:

```bash
python -m build
pip install dist/relatipy-0.0.0-py3-none-any.whl
```

To publish to PyPI use `twine`:

```bash
python -m build
twine upload dist/*
```
