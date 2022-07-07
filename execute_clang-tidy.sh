#!/bin /bash

clang-tidy -header-filter=.* -extra-arg=-std=c++20 src/* > errors.log