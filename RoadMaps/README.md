#Road Maps
To combine road maps without removing duplicates use:
```
awk 'FNR==1{print ""}1' *.txt > finalfile.txt
```