import os
import numpy as np
import ctypes
from numpy.ctypeslib import ndpointer


def FwdBwdAlg_cpp(initPi, transPi, SoftEv, order='C'):

    if not hasEigenLibReady:
        raise ValueError("Cannot find library %s. Please recompile."
                         % (libfilename))
    if order != 'C':
        raise NotImplementedError("LibFwdBwd only supports row-major order.")

    T, K = SoftEv.shape

    # Prep inputs
    initPi = np.asarray(initPi, order=order)
    transPi = np.asarray(transPi, order=order)
    SoftEv = np.asarray(SoftEv, order=order)

    # Allocate outputs
    resp = np.zeros((T, K), order=order)
    resp_pair = np.zeros((K, K), order=order)
    marg_pr_seq = np.zeros((1, 1), order=order)

    # Execute C++ code (fills in outputs in-place)
    lib.FwdBwdAlg(initPi, transPi, SoftEv, resp, resp_pair, marg_pr_seq, K, T)

    return resp, resp_pair, marg_pr_seq


path = os.getcwd()
libpath = path + '/ChromA/util/'
libfilename = 'libfwdbwdcpp.so'
hasEigenLibReady = True

try:
    lib = ctypes.cdll.LoadLibrary(os.path.join(libpath, libfilename))
    lib.FwdBwdAlg.restype = None
    lib.FwdBwdAlg.argtypes = \
        [ndpointer(ctypes.c_double),
         ndpointer(ctypes.c_double),
         ndpointer(ctypes.c_double),
         ndpointer(ctypes.c_double),
         ndpointer(ctypes.c_double),
         ndpointer(ctypes.c_double),
         ctypes.c_int, ctypes.c_int]

except OSError:
    # No compiled C++ library exists
    print("Failed to Load Cpp Core")
    hasEigenLibReady = False
