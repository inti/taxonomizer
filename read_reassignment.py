import numpy as np
import scipy
import h5py
import ctime

# Make class to update model

class ReadAssignment:
	def __init__(self,hdf5_file,niter = 100):
		# init parameters of the model
		self.hdf5_file = hdf5_file
		self.nR,self.nG = self.hdf5_file["Q_all"].shape
		# initialize tables for parameters on disk
		self.params = dict()
	#def initialize_parameters():
                if "pathoscope2_parameters" in self.hdf5_file:
                        del self.hdf5_file["pathoscore2_parameters"]
		
                pathoscore2_parameters = self.hdf5_file.require_group("pathoscore2_parameters")
                self.params["pi"] = pathoscore2_parameters.create_dataset("pi",(self.nG,),dtype='float')
                self.params["theta"] = pathoscore2_parameters.create_dataset("theta",(self.nG,),dtype='float')
                self.params["a"] = pathoscore2_parameters.create_dataset("a",(self.nG,),dtype='float')
                self.params["b"] = pathoscore2_parameters.create_dataset("b",(self.nG,),dtype='float')
                self.params["delta"] = pathoscore2_parameters.create_dataset("delta",(self.nR,self.nG),dtype='float')


h5 = h5py.File("mytestfile.hdf5", "r+")

ra = ReadAssignment(h5,niter=50)


#  pathoscope2

# initialize parameters

ra.params["pi"][:] = 1/float(ra.nG)
ra.params["theta"][:] = 1/float(ra.nG)
ra.params["a"][:] = 1
ra.params["b"][:] = 1


print ctime(),"Done"
