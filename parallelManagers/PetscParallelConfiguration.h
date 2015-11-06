#ifndef _PETSC_PARALLEL_CONFIGURATION_H_
#define _PETSC_PARALLEL_CONFIGURATION_H_

#include "../Parameters.h"
#include "../Definitions.h"
#include <mpi.h>


/** Class used to set parameters relevant to the parallel distribution. All functions modify the
 * parameters instance
 */
class PetscParallelConfiguration {

    private:

        Parameters & _parameters;   //! Reference to the parameters

        /** Locates the six neighbors of the current process
         */
        void locateNeighbors();

        /** Computes the indices of the current subdomain and stores the data in the parameters
         */
        void createIndices();

        /** Returns the rank of the process with the indices provided
         * @param i Intex in the X directon
         * @param j Intex in the Y directon
         * @param k Intex in the Z directon
         * @return Rank of the process with that index, assuming that they are ordered
         * lexicographically, or MPI_PROC_NULL if outside the domain
         */
        int computeRankFromIndices(int i, int j, int k) const;

        /** Compute local sizes and sizes in all directions. Requires deallocation of sizes
         */
        void computeSizes();

        /** Deletes the arrays allocated in the parameters. To be called in the destructor of this
         * class
         */
        void freeSizes();

    public:

        /** Constructor
         * @param parameters Reference to the parameters
         */
        PetscParallelConfiguration(Parameters & parameters);

        /** Destructor */
        ~PetscParallelConfiguration();
};

#endif
