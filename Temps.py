bmin = 0.0125
bmax = 1.37
N = 61
b = [(bmax-bmin)*i/N + bmin for i in range(N+1)]

b = list(set(b))
b = sorted(b,key=lambda x:x)

fout = open('temps.dat','w')
for i in b:
    fout.write('%0.6f\n' % i)
fout.close()
