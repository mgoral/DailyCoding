#!/bin/bash

# usage:
# build <quiz directory name>

SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
QUIZ_DIR="$1"

pushd "${SCRIPT_DIR}/${QUIZ_DIR}"
g++ main.cpp -std=c++14 && ./a.out
popd
