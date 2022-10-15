#!/bin/bash

case "$1" in
    -a | --all ) bazelisk run networks -- "$2"; python sandbox.py all; feh -. graphics/test_graphic.png &;;
    -b | --build ) bazelisk build networks;;
    -r | --run ) bazelisk run networks  -- "$2";;
    -d | --draw ) python sandbox.py all; feh -. graphics/test_graphic.png;;
    -s | --show ) feh -. graphics/test_graphic.png &;;
    -sf | --show_file ) feh -. graphics/"$2".png &;;
    -g | --git ) git add .; git commit -m "$2"; git push origin main;;
    -u | --upload ) scp -r ./graphics/test_graphic.png senghaas@LIN5G.tphys.uni-heidelberg.de:~/WWW/research/files/network.png;;
    -v | --venv ) source .venv/bin/activate;;
esac
