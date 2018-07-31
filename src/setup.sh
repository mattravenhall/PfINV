#!/bin/bash

# Check for, then install, nucmer
if ! [ -x "$(command -v nucmer)" ]; then sudo apt-get install mummer; fi
if ! [ -x "$(command -v show-coords)" ]; then sudo apt-get install mummer; fi