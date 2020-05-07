#!/usr/bin/awk -f
BEGIN{sum=0; sumsqd=0; n=0}
{sum+=$1; sumsqd+=$1*$1; n++}
END{print sum/n" "sqrt((sumsqd/n-(sum/n)*(sum/n))*(n/(n-1)))}
