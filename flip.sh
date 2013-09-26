#!/usr/bin/env bash

echo "dof  u jacobi1     chebyshev1  ssor1       jacobi4     chebyshev4  ssor4       jacobi16    chebyshev16 ssor16" > $1.bak
psc -r < $1 | sc -W% - >> $1.bak
mv $1.bak $1



