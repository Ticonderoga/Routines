#!/bin/sh

ipython -pylab "DiagHumid.py"

echo "

------------------
(program exited with code: $?)" 		


unlink $0
