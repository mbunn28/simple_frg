import matplotlib.pyplot as plt

def readFile(fileName):
        fileObj = open(fileName, "r") #opens the file in read mode
        words = fileObj.read().splitlines() #puts the file into an array
        fileObj.close()
        return words

def main():
    lambda_vals = readFile("lambda_vals")
    tp_vals = readFile("tp_vals")

    plt.plot(tp_vals,lambda_vals)
    plt.savefig("fig.png")
    return  

main()