import numpy as np
import matplotlib.pyplot as plt
import glob

def main():
    fnames = glob.glob("obs_E_*")
    fnames = sorted(fnames,key=lambda x:float(x.split("_")[2]))
    plotme = []
    for f in fnames:
        d = np.loadtxt(f)
        plotme.append([float(f.split("_")[2]),d.mean(),d.std()/(len(d)**0.5)])
    plotme = np.array(plotme)
    plt.errorbar(1.0/plotme[:,0],plotme[:,1],plotme[:,2])
    plt.show()

def main2():
    names = ["L6.dat","L8.dat","L10.dat"] 
    for f in names:
        s = f.split('L')[1]
        s = s.split('.')[0]
        s = int(s)
        d = np.loadtxt(f)
        plt.plot(1./d[:,0],d[:,2]/(s*s),label="%d"%s)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main2()
