from matplotlib.pyplot import *

def genDiagram(fileName, x, xlab, ys, ylab, yleg):
    figure(dpi=300)
    for i in range(len(ys)):
      plot(x, ys[i], label=yleg[i])
    xlabel(xlab)
    ylabel(ylab)
    title('Speciation Diagram')
    legend()
    tight_layout()
    savefig(fileName)