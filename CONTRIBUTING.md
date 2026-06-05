# Contributing to amber-codon-scanner

Thank you for your interest in contributing!

## Reporting bugs

Please open an issue on GitHub with:
- A description of the bug
- A minimal reproducible example (DNA sequence and expected output)
- Your Python version and operating system

## Suggesting features

Open an issue describing:
- The feature you would like
- The biological use case it addresses
- Any relevant literature on pyrrolysine or amber codon biology

## Contributing code

1. Fork the repository
2. Create a new branch: `git checkout -b my-feature`
3. Make your changes
4. Run the tests: `pytest`
5. Push your branch and open a pull request

## Development setup

```bash
git clone https://github.com/CameronCat/amber-codon-scanner
cd amber-codon-scanner
pip install -e ".[dev]"
pytest
```

## Code style

- Follow PEP 8
- Add docstrings to new functions
- Add tests for new functionality
- Keep biological claims grounded in published literature with citations
- Be explicit about heuristic vs validated methods

## Questions

Open an issue on GitHub.