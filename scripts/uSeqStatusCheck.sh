#!/bin/bash

echo "Running:"
echo "$(qstat | awk '$5 ~ /^r$/ {print $3} ' | sort | uniq -c)"

echo "Waiting:"
echo "$(qstat | awk '$5 ~ /^qw$/ {print $3} ' | sort | uniq -c)"

echo "Holding:"
echo "$(qstat | awk '$5 ~ /^hqw$/ {print $3} ' | sort | uniq -c)"
