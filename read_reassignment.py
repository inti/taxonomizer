import numpy as np
import scipy
import h5py
import ctime

h5 = h5py.File("mytestfile.hdf5", "r+")

Q_all_row_sums = h5.create_dataset("Q_all_row_sums",(h5["Q_all"].shape[0],),dtype="float")

Q_all_row_sums[:] = np.sum(h5["Q_all"],1)

# make the x latent variable : X = h5["Q_all"] / Q_all_row_sums[()][:,np.newaxis]








print ctime(),"Done"
