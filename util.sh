#!/bin/bash

case "$1" in
    -a | --all ) bazel run networks -- "$2"; python sandbox.py all; feh -. graphics/test_graphic.png &;;
    -b | --build ) bazel build networks;;
    -r | --run ) bazel run networks  -- "$2";;
    -d | --draw ) python sandbox.py all;;
    -s | --show ) python sandbox.py all; feh -. graphics/test_graphic.png &;;
    -sf | --show_file ) feh -. graphics/"$2".png &;;
    -g | --git ) git add .; git commit -m "$2"; git push origin main;;
esac
