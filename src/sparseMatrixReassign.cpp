/* sparseMatrixReassign. This MEX function efficiently updates values in a
 * sparse matrix. This circumvents the need to re-allocate memory and
 * re-initialize memory every time matrix entries need to be updated while
 * advancing in time.
*/


extern "C"
{
    #include "mex.h"
    #include "matrix.h"
}

void sparseMatrixReassignFunction(mxDouble* sparseMatrixVals, const mxDouble* newMatrixVals, const mxInt32* cooToCscIndices, const mwIndex nnz) {
  
    for(mwIndex j=0; j<nnz; j++) {
    sparseMatrixVals[cooToCscIndices[j]-1] = newMatrixVals[j];
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nlhs > 0 || nrhs != 3) return;
    if(!mxIsSparse(prhs[0])) return;

    mxArray* nonConstSparseMatrix = const_cast<mxArray*>(prhs[0]);

    mwSize numColumns;
    mwIndex *columnIndexVector;
    mwIndex nnz;
    mxDouble *sparseMatrixVals, *newMatrixVals;

    numColumns = mxGetN(nonConstSparseMatrix);
    columnIndexVector = mxGetJc(nonConstSparseMatrix);
    sparseMatrixVals = mxGetDoubles(nonConstSparseMatrix);
    nnz = columnIndexVector[numColumns];

    newMatrixVals = mxGetDoubles(prhs[1]);

    mxInt32* COO2CSCIndices = mxGetInt32s(prhs[2]);

    sparseMatrixReassignFunction(sparseMatrixVals, newMatrixVals, COO2CSCIndices, nnz); 

}