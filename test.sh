#!/bin/bash

# Find all *-test.c files

test_files=$(ls ./src/*-test.c)

# Compile each test file
for test_file in $test_files;
do
    gcc -o ${test_file%.c}.out $test_file

    if [ $? -eq 0 ];
    then
        echo "testing ${test_file%.c}"
        ./${test_file%.c}.out
    else
        echo "Compilation failed for $test_file"
    fi
    rm ./${test_file%.c}.out
done
echo "Done"
