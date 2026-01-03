#!/bin/bash

# Run the build command in the venv
python setup.py build_ext --inplace

python setup.py clean --all