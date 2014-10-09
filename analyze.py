# Analysis script for 2D Ising Model

import tarfile
import numpy as np
import sys

def get_temps(s1):
    betas = []
    for i in s1.getnames():
        try:
            t = i.split('obs_ratio_')[1]
            if t not in betas:
                betas.append(t)
        except:
            continue
    return betas

def good_reps(data):
    beta = get_temps(data)
    fmask = 'L%s/R%s/r%03d/obs_ratio_%s'
    f_reps = {}
    for names in data.getnames():
        L = names.split('/')[0].split('L')[1]
        if L not in f_reps:
            f_reps[L] = {}
        R = names.split('/')[1].split('R')[1]
        if R not in f_reps[L]:
            f_reps[L][R] = []
        r = names.split('/')[2].split('r')[1]
        if r not in f_reps[L][R]:
            f_reps[L][R].append(r)
    names = data.getnames()
    # They are built in steps of one, each corresponding to a row of the LxL system
    for L in f_reps.keys():
        for R in f_reps[L].keys():
            for B in beta:
                for zz in range(int(L)):
                    if fmask % (L,R,zz,B) not in names:
                        try:
                            del f_reps[L][R]
                            break
                        except KeyError:
                            break
        if len(f_reps[L])==0:
            del f_reps[L]
    return f_reps

def get_union(s1,s2):
    for L in s1.keys():
        for R in s1[L].keys():
            if R not in s2[L]:
                del s1[L][R]
    for L in s2.keys():
        for R in s2[L].keys():
            if R not in s1[L]:
                del s2[L][R]
    return s1

def dump_ent(L,B,l1,l2):
    # Both are two dimensiona arrays, 1st index is disorder, 2nd index is region size
    # Working at a particular L and B, which are strings
    fout = open('L%sB%s' % (L,B),'w')
    for i in l1:
        for j in i:
            fout.write('%0.8f ' % j)
        fout.write('\n')
    fout.write('\n')
    for i in l2:
        for j in i:
            fout.write('%0.8f ' % j)
        fout.write('\n')
    fout.close()
    fout = open('MI_L%sB%s' % (L,B),'w')
    for i1,i in enumerate(l1):
        for j1,j in enumerate(i):
            val = l1[i1][j1] + l2[i1][len(l2[i1])-1-j1] - l1[i1][-1]
            val += l2[i1][j1] + l1[i1][len(l2[i1])-1-j1] - l2[i1][-1]
            val /= 2.0
            fout.write('%0.8f ' % val)
        fout.write('\n')
    fout.write('\n')
    fout.close()

def logmean(d1,d2):
    # Average the files d1 and d2
    # d1 is the summed data, where d2 is the multiplier for each block
    maxval = d2.max()
    t = 0
    n = len(d1)
    for i in range(len(d1)):
        t += d1[i] * np.exp(d2[i]-maxval)
    return np.log(t/n) + maxval

def main():
    data_f = tarfile.open('2014-10-06-CleanTest/data.tar','r')
    data_r = tarfile.open('2014-10-06-CleanTest/data.tar','r')
    # print data_f.getnames() # Of the form L###/R###/r###/obs_ratio_#.######
    # L - size, R - replica, r - ratio, last is beta
    # First, we have to find all the sizes we have to check
    # For each size, check what replicas we have to check
    # This is the set of good sizes and replicas for which we have data
    betas = get_temps(data_f)
    fmask = 'L%03d/r%03d/obs_ratio_%s'
    facmask = 'L%03d/r%03d/fac_%s'
    all_data = {}
    mid_vals = {}
    minset = [6,8,10]
    for L in minset:
        for B in betas:
            # These will contain the entropies, we will only do error as disorder error
            # With both of these we can also calculate mutual information
            t_ent_f = [0 for i in range(2*int(L)+1)]
            t_ent_r = [0 for i in range(2*int(L)+1)]
            for zz in range(2*int(L)):
                df = np.loadtxt(data_f.extractfile(fmask % (L,zz,B)))
                dffac = np.loadtxt(data_f.extractfile(facmask % (L,zz,B)))
                t_ent_f[zz+1] = t_ent_f[zz] - logmean(df,dffac)
                dr = np.loadtxt(data_r.extractfile(fmask % (L,zz,B)))
                drfac = np.loadtxt(data_r.extractfile(facmask % (L,zz,B)))
                t_ent_r[zz+1] = t_ent_r[zz] - logmean(dr,drfac)
            #dump_ent(L,B,t_ent_f,t_ent_r)
            ent = [0 for i in range(2*int(L)+1)]
            #ent2 = [0 for i in range(int(L)+1)]
            mi = [0 for i in range(2*int(L)+1)]
            #mi2 = [0 for i in range(int(L)+1)]
            for zz in range(2*int(L)):
                ent[zz+1] = t_ent_f[zz+1]
                #ent2[zz+1] = t_ent_f[zz+1]**2
                mi[zz+1] = (t_ent_f[zz+1] + t_ent_r[2*int(L)-zz-1] - t_ent_f[2*int(L)] + t_ent_r[zz+1] + t_ent_f[2*int(L)-zz-1] - t_ent_r[2*int(L)])/2.0
                #mi2[zz+1] = ((t_ent_f[zz+1] + t_ent_r[int(L)-zz-1] - t_ent_f[int(L)] + t_ent_r[zz+1] + t_ent_f[int(L)-zz-1] - t_ent_r[int(L)])/2.0)**2.0
            ent = np.array(ent)
            #ent2 = np.array(ent2)
            mi = np.array(mi)
            #mi2 = np.array(mi2)
            #all_data[(L,B)] = [ent,((ent2-ent**2)/N)**0.5,mi,((mi2-mi**2)/N)**0.5]
            all_data[(L,B)] = [ent,mi]
            #mid_vals[(L,B)] = [ent[int(L)/2],((ent2[int(L)/2]-ent[int(L)/2]**2)/N)**0.5,mi[int(L)/2],((mi2[int(L)/2]-mi[int(L)/2]**2)/N)**0.5]
            mid_vals[(L,B)] = [ent[int(L)],mi[int(L)]]
    pa = sorted(all_data,key=lambda x:x[1])
    sorted(pa,key=lambda x:x[0])
    for i in pa:
        fout = open('L%s.dat' % i[0],'a')
        #fout.write('%s %0.8f %0.8f %0.8f %0.8f\n' % (i[1],mid_vals[i][0],mid_vals[i][1],mid_vals[i][2],mid_vals[i][3]))
        fout.write('%s %0.8f %0.8f\n' % (i[1],mid_vals[i][0],mid_vals[i][1]))
        fout.close()

if __name__ == "__main__":
    main()
