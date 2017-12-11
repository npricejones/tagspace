# Write a function to grab all instances that use a particular set of generation functions
# Write a function to grab all instances within a certain time period
# Write a function to check the attributes of a given group or dataset against a provided filter

def external_validation(labels_true,labels_pred,metric):
    Nratio = labels_true.shape[0]
    Niters = labels_true.shape[1]
    success = np.zeros((Nratio,Niters))
    for r in range(Nratio):
        for n in range(Niters):
            success[r][n] = metric(labels_true[r][n],labels_pred[r][n])
            
    return success
            