#!/bin/bash
cd /data
# Check if the CDK_JAR environment variable is set
if [ -z "$CDK_JAR" ]; then
    echo "Error: CDK_JAR environment variable is not set."
    echo "Please set the CDK_JAR environment variable to the path of the CDK library jar file."
    exit 1
fi

# Check if the CDK library jar file exists
if [ ! -f "$CDK_JAR" ]; then
    echo "Error: The CDK library jar file specified by CDK_JAR does not exist."
    echo "Path specified: $CDK_JAR"
    exit 1
fi

# Check if the Java source file exists
if [ ! -f "CDKTest.java" ]; then
    echo "Error: CDKTest.java does not exist in the current directory."
    exit 1
fi

# Compile the Java program
echo "Compiling CDKTest.java..."
javac -cp .:"$CDK_JAR" CDKTest.java
if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

# Run the Java program
echo "Running CDKTest..."
java -cp .:"$CDK_JAR" CDKTest
