import numpy as np
import sys
import copy


class ParamStorage(object):
    """
    Minimal Parameter Storage Unit with some Checks
    Inspired by: Michael Hughes's ParamBag
    """

    def __init__(self, K=0, doCollapseK1=False, **kwargs):
        self.K = K
        self.D = 0
        if sys.version_info > (3, 0):
            for key, val in kwargs.items():
                setattr(self, key, val)
        else:
            for key, val in kwargs.iteritems():
                setattr(self, key, val)
        self._FieldDims = dict()
        self.doCollapseK1 = doCollapseK1

    def setField(self, key, rawArray, dims=None):
        # Parse dims tuple
        if dims is None and key in self._FieldDims:
            dims = self._FieldDims[key]
        else:
            self._FieldDims[key] = dims
        # Parse value as numpy array
        setattr(self, key, self.parseArr(rawArray, dims=dims, key=key))

    def parseArr(self, arr, dims=None, key=None):
        K = self.K
        D = self.D
        arr = np.asarray(arr, dtype=np.float64)
        # Verify shape is acceptable given expected dimensions
        if dims is not None and isinstance(dims, str):
            dims = (dims)  # force to tuple
        expectedShape = self._getExpectedShape(dims=dims)
        if self.doCollapseK1:
            if arr.shape not in self._getAllowedShapes(expectedShape):
                self._raiseDimError(dims, arr, key)
            # Squeeze into most economical shape possible
            #  e.g. (3,1) --> (3,),  (1,1) --> ()
            if K == 1 or D == 1:
                arr = np.squeeze(arr)
        else:
            if arr.shape != expectedShape:
                self._raiseDimError(dims, arr, key)
        if arr.ndim == 0:
            arr = np.float64(arr)
        return arr

    def _getExpectedShape(self, key=None, dims=None):
        ''' Returns tuple of expected shape, given named dimensions.

        Example
        -------
        >>> PB = ParamStorage(K=3, D=2)
        >>> PB._getExpectedShape(dims=('K','K'))
        (3, 3)
        >>> PB._getExpectedShape(dims=('K','D','D'))
        (3, 2, 2)
        '''
        if key is not None:
            dims = self._FieldDims[key]
        if dims is None:
            expectShape = ()
        else:
            shapeList = list()
            for dim in dims:
                if isinstance(dim, int):
                    shapeList.append(dim)
                else:
                    shapeList.append(getattr(self, dim))
            expectShape = tuple(shapeList)
        return expectShape
