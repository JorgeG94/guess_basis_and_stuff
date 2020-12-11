#!/bin/bash

#export JSON_PATH=${HOME}/programs/prep/json/single_include
export JSON_PATH=${HOME}/Programs/Install/nlohmann-json
g++ -std=c++11 -I${JSON_PATH}/include json.cpp
