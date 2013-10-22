/*
 * FEL_incremental_gpu.h
 *
 *  Created on: Sep 15, 2012
 *      Author: abon
 */

#ifndef FEL_INCREMENTAL_GPU_H_
#define FEL_INCREMENTAL_GPU_H_


__device__ float GetElement(const Matrix A, int row, int col)
{
    return A.elements[row * A.stride + col];
}
// Set a matrix element
__device__ void SetElement(Matrix A, int row, int col,
                           float value)
{
    A.elements[row * A.stride + col] = value;
}



#endif /* FEL_INCREMENTAL_GPU_H_ */
