#!/bin/bash

cd ./src/

javac ./*.java

cd ..

mv ./src/*.class ./bin/
mv ./src/org/knowceans/gibbstest/*.class ./bin/org/knowceans/gibbstest/
