

import sys


def fib(n):
    if n < 0: return -1
    if n== 0: return 0
    if n== 1: return 1
    return fib(n-1)+fib(n-2)

n=int(sys.argv[1])
print(f"n={n}")
for i in range(n+1):

    print(f"fib({i})={fib(i)}")
