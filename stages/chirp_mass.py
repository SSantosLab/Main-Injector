import pickle
import numpy as np

def get_chirp_mass(distmean, diststd, area90, area50):
    test_event = [[distmean, diststd, area90, area50]] # distmean, diststd, area90, area50

    with open('chirp_mass.dat', 'rb') as f:
        clf = pickle.load(f)

    pdf = clf.predict_proba(test_event)[0]
    if clf.missing_bins.shape[0]:
        for c in clf.missing_bins:
            pdf = np.insert(pdf, c-1, 0.)
    expv = np.sum(pdf*clf.midpoints)
    std = np.sqrt(np.sum(pdf*clf.midpoints*clf.midpoints) - expv**2)

    return [round(expv,1), round(std,1)]