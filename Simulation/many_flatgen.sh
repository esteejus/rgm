#!/bin/bash

for i in {1..100}
do
  root -b -q "generate_neutrons.C(\"lundfiles/flatgen_$i.txt\", 100000, 35, 135 )"
done
