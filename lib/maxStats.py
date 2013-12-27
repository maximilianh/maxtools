class Stats:
    """ just a bag of values, to group together results """
    pass

    #def __repr__(self):
        #lines = []
        #lines.append("TP   %f" % self.TP)
        #lines.append("FP   %f"  % self.FP)
        #lines.append("TN   %f"  % self.TN)
        #lines.append("FN   %f"  % self.FN)
        #lines.append("Sens %f"  % self.Sens)
        #lines.append("Spec %f"  % self.Spec)
        #lines.append("PPV  %f"  % self.PPV)
        #lines.append("CC   %f"  % self.CC)
        #return "\n".join(lines)

def divide(top, bottom):
    if bottom!=0:
        return float(top) / float(bottom)
    else:
        return 100000000000000

def genomeLevelCorrelation(genomeSize, predictionSize, referenceSize, overlap):
    """ calculate correlation based on overlap of a prediction and a reference genome track """
    TP = overlap
    FP = sizeA - overlap
    FN = sizeB - overlap
    TN = genomeSize - 
    
def calcCorrelation(TP, FP, FN, TN):
    """ pearson correlation """
    stats.TP = TP
    stats.FP = FP
    stats.TN = TN
    stats.FN = FN

    stats.Sens = divide(TP , (TP + FN))
    stats.Spec = divide(TN , (TN + FP))
    stats.PPV  = divide(TP , (TP + FP))      # PRECISION aka Positive Predictive Value
    # Precision measures the proportion of the claimed true functional sites that are indeed true functional sites.
    stats.PC   = divide(TP , (TP + FP) )
    # Accuracy  measures the proportion of predictions, both for true functional sites and false functional sites that are correct. 
    stats.Acc  = divide((TP + TN) , (TP + FP + FN + TN))

    CC_top = TP * TN - FN * FP
    CC_bottom = math.sqrt((TP+FN)*(TN+FP)*(TP+FP)*(TN+FN))
    stats.CC = divide( CC_top , CC_bottom )

