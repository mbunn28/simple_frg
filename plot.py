import matplotlib.pyplot as plt
import numpy as np

def readFile(fileName):
        fileObj = open(fileName, "r") #opens the file in read mode
        words = fileObj.read().splitlines() #puts the file into an array
        fileObj.close()
        return words

def main():
    lambda_vals_phases = readFile("lambda_vals")
    tp_vals = readFile("tp_vals")
    tp_vals = np.array(tp_vals,dtype=float)
    tp_vals = -np.round(tp_vals, 2)
    lambda_vals = lambda_vals_phases[0::2]
    lambda_vals = np.array(lambda_vals,dtype=float)
    phases = lambda_vals_phases[1::2]
    phases = np.array(phases)

    sc = phases == 'SC'
    afm = phases == 'AFM'
    fm = phases == 'FM'

    # print(lambda_vals)
    # print(phases)
    # print(tp_vals)
    plt.semilogy(tp_vals[sc],lambda_vals[sc],'kx')
    plt.semilogy(tp_vals[afm],lambda_vals[afm],'bx')
    plt.semilogy(tp_vals[fm],lambda_vals[fm],'rx')
    plt.xlabel(r"$-t'$ / $t$")
    plt.ylabel(r'$\Lambda_C$')
    plt.legend(["SC","AFM","FM"])
    plt.savefig("fig.png")
    return  

main()