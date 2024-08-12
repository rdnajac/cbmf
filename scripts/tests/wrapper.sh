#!/bin/bash
#
## Wrapper script to test sourcing another script.

# source dontrunme.sh

# ./dontrunme.sh is not the same as the next line!!!
. dontrunme.sh

# but this this same
# . ./dontrunme.sh


SUCCESS
