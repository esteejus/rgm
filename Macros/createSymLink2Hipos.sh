#!/bin/bash
SOURCE=$1
DIR=$2
find $SOURCE -name \*.hipo -exec ln -vs "{}" $DIR ';'
