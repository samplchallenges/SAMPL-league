#!/bin/bash -e

if id "webapp" &>/dev/null; then
    :
else
    sudo useradd -m webapp
fi
