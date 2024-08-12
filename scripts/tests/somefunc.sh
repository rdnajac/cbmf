#!/bin/bash
#
## Just some function

myfunc() 
{
    echo "Hello from myfunc"
}

# prints message:
# ./somefunc.sh && myfunc
#
# this works the same as using 'source' and then calling the function
# but 'source' is not POSIX compliant and '.' does the same thing

# does not print message:
# ./somefunc.sh
# myfunc
# ./somefunc.sh myfunc

# if uncomment the following line
# "$@"

# then we can call the function with:
# ./somefunc.sh myfunc

# the "$@" tells the script to expect positional parameters
# and to pass them to the script when it is called
#
# basically it lets the script know to process arguments when it is called
# (I think)
