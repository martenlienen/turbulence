#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_

#include <string>
#include <petscksp.h>
#include "Definitions.h"
#include "Parameters.h"

class Configuration{
    private:
        std::string _filename;
        int _dim;

    public:
        Configuration();
        Configuration(const std::string & filename);
        void setFileName (const std::string & filename);
        void loadParameters(Parameters & parameters,
                            const MPI_Comm & communicator = PETSC_COMM_WORLD);
};

#endif
