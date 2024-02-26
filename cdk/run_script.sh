#!/bin/bash

# Compile the Java program
javac -cp .:$CDK_JAR CDKTest.java

# Run the Java program
java -cp .:$CDK_JAR CDKTest
