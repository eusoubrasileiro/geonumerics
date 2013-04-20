#!/bin/bash
# convert to little endian all intel based machines
# suswapbytes < data.su format=1 > data.su.swapped

# not easy to find in su documentation, easily convert from su to segy
segywrite < $1 tape=out.sgy

# the 'binary' 400 bytes and 'header' 3200 text files are created by segywrite


