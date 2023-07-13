#!/bin/bash

readonly EXEC_DIR=$(pwd)

cd "bayes-mapper.pyのあるディレクトリ"

readonly SCRIPT_DIR=$(pwd)

# enter venv
. .venv/bin/activate

# execute AU-aLRT
cd $EXEC_DIR
python3 "${SCRIPT_DIR}/bayes-mapper.py" $@
# exit venv
deactivate
