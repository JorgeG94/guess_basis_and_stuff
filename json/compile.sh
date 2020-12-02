#!/bin/bash

#export JSON_PATH=${HOME}/programs/prep/json/single_include
export JSON_PATH=${HOME}/programs/install/json/include
g++ -std=c++11 -L${JSON_PATH} json.cpp
