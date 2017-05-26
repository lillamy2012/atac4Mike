import sys
from collections import Counter


cnt = Counter()
for line in sys.stdin:
    lines = line.split()
    cnt[abs(int(lines[0]))] +=1
    # sorts items
items = sorted(cnt.items(), key=lambda x: int(x[0]))

for k, repetitions in items:
    print k,'\t', repetitions
