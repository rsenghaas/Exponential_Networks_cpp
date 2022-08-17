#!/bin/bash

case "$1" in
    -a | --all ) bazel run networks -- "$2"; python sandbox.py all;;
    -b | --build ) bazel build networks;;
    -r | --run ) bazel run networks  -- "$2";;
    -d | --draw ) python sandbox.py all;;
    -s | --show ) python sandbox.py all; feh -. graphics/test_graphics.png &;;
    -sf | --show_file ) feh -. graphics/"$2".png &;;
esac
