#ifndef MOLDB_CUDA_H
#define MOLDB_CUDA_H

#include "util.h"

#include <cuda_runtime.h>

namespace MolDB {

// This will output the proper CUDA error strings in the event that a CUDA host call returns an error
#define checkCudaErrors(err)  __checkCudaErrors (err, __FILE__, __LINE__)

  inline std::string __checkCudaErrors(cudaError err, const char *file, const int line )
  {
    if(cudaSuccess != err)
      return make_string(file, "(", line, "): CUDA Runtime API error ", static_cast<int>(err), ": ", cudaGetErrorString(err), ".");
    return std::string();
  }

// This will output the proper error string when calling cudaGetLastError
#define getLastCudaError(msg)      __getLastCudaError (msg, __FILE__, __LINE__)

  inline std::string __getLastCudaError(const char *errorMessage, const char *file, const int line )
  {
    cudaError_t err = cudaGetLastError();
    if (cudaSuccess != err)
      return make_string(file, "(", line, "): getLastCudaError() CUDA error: ", errorMessage, ": (", static_cast<int>(err), ") ", cudaGetErrorString(err), ".");
    return std::string();
  }

}
#endif
