
import os
import matplotlib.pyplot as plt

# Find out the name of the directory containing this file
__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

datafiles = [fname for fname in os.listdir(__location__) if "data" in fname]

for fname in datafiles:
    # Generate x and y's
    with open(os.path.join(__location__, fname)) as f:
        lines = f.readlines()
        # First line gives ranges
        ranges = lines[0].split()
        x = [float(line.split()[0]) for line in lines[1:]]
        y = [float(line.split()[1]) for line in lines[1:]]
        
        plt.scatter(x, y, marker = "+", c = "tab:blue")

        if len(ranges) > 2:
            xmin = float(ranges[0])
            xmax = float(ranges[1])
            ymin = float(ranges[2])
            ymax = float(ranges[3])
            plt.xlim([xmin, xmax])
            plt.ylim([ymin, ymax])
            
        xlabel = ranges[-2]
        ylabel = ranges[-1]
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        plotname = os.path.join(__location__, "plot-{}.png".format(fname[5:]))
        plt.savefig(plotname)
        plt.clf() # clear
        
