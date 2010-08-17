#!/bin/sh
dataset=$1
shift
timelimit="Inf"
while test -n "${1}"; do
    case ${1} in
        --time-limit)
            timelimit=${2}
            shift
            break;;
    esac
    shift
done
R --vanilla <<EOF
library(retistruct)
retistruct.cli("$dataset", $timelimit)
EOF

