#include <float.h> // To get the max double
#include <math.h>
#include "SORSolver.h"

SORSolver::SORSolver(FlowField & flowField, const Parameters & parameters):
    LinearSolver(flowField, parameters){}

void SORSolver::solve() {
    // Solves the heat equation
    FLOAT resnorm = DBL_MAX, tol = 1e-4;

    double omg = 1.7;
    int iterations = -1;
    int it = 0;

    int nx = _flowField.getNx(), ny = _flowField.getNy(), nz = _flowField.getNz();
    ScalarField & P = _flowField.getPressure();
    if (_parameters.geometry.dim == 3){
        do {
            for (int k = 2; k < nz + 2; k++){
                for (int j = 2; j < ny + 2; j++){
                    for (int i = 2; i < nx + 2; i++){
                        const FLOAT dx_0 = _parameters.meshsize->getDx(i  ,j  ,k  );
                        const FLOAT dx_M1= _parameters.meshsize->getDx(i-1,j  ,k  );
                        const FLOAT dx_P1= _parameters.meshsize->getDx(i+1,j  ,k  );
                        const FLOAT dy_0 = _parameters.meshsize->getDy(i  ,j  ,k  );
                        const FLOAT dy_M1= _parameters.meshsize->getDy(i  ,j-1,k  );
                        const FLOAT dy_P1= _parameters.meshsize->getDy(i  ,j+1,k  );
                        const FLOAT dz_0 = _parameters.meshsize->getDz(i  ,j  ,k  );
                        const FLOAT dz_M1= _parameters.meshsize->getDz(i  ,j  ,k-1);
                        const FLOAT dz_P1= _parameters.meshsize->getDz(i  ,j  ,k+1);

                        const FLOAT dx_W = 0.5*(dx_0+dx_M1);
                        const FLOAT dx_E = 0.5*(dx_0+dx_P1);
                        const FLOAT dx_S = 0.5*(dy_0+dy_M1);
                        const FLOAT dx_N = 0.5*(dy_0+dy_P1);
                        const FLOAT dx_B = 0.5*(dz_0+dz_M1);
                        const FLOAT dx_T = 0.5*(dz_0+dz_P1);

                        const FLOAT a_W  =  2.0/(dx_W*(dx_W+dx_E));
                        const FLOAT a_E  =  2.0/(dx_E*(dx_W+dx_E));
                        const FLOAT a_N  =  2.0/(dx_N*(dx_N+dx_S));
                        const FLOAT a_S  =  2.0/(dx_S*(dx_N+dx_S));
                        const FLOAT a_T  =  2.0/(dx_T*(dx_T+dx_B));
                        const FLOAT a_B  =  2.0/(dx_B*(dx_T+dx_B));
                        const FLOAT a_C  = -2.0/(dx_E*dx_W)-2.0/(dx_N*dx_S)-2.0/(dx_B*dx_T);

                        P.getScalar(i,j,k) = omg/a_C*(
                                               _flowField.getRHS().getScalar(i,j,k)
                                               - a_W*P.getScalar(i-1,j,k) - a_E*P.getScalar(i+1,j,k)
                                               - a_S*P.getScalar(i,j-1,k) - a_N*P.getScalar(i,j+1,k)
                                               - a_B*P.getScalar(i,j,k-1) - a_T*P.getScalar(i,j,k+1) )
                                            + (1.0-omg) * P.getScalar(i,j,k);
                    }
                }
            }

            for ( int j = 2; j < ny + 2; j++ ) {
                for (int k = 2; k < nz + 2; k++ ) {
                    P.getScalar(1,j,k) = P.getScalar(2,j,k);
                    P.getScalar(nx+2,j,k) = P.getScalar(nx+1,j,k);
                }
            }

            for ( int i = 2; i < nx + 2; i++ ) {
                for (int k = 2; k < nz + 2; k++ ) {
                    P.getScalar(i,1,k) = P.getScalar(i,2,k);
                    P.getScalar(i,ny+2,k) = P.getScalar(i,ny+1,k);
                }
            }

            for ( int i = 2; i < nx + 2; i++ ) {
                for (int j = 2; j < ny + 2; j++) {
                    P.getScalar(i,j,1) = P.getScalar(i,j,2);
                    P.getScalar(i,j,nz+2) = P.getScalar(i,j,nz+1);
                }
            }

            resnorm = 0;
            for (int k = 2; k < nz + 2; k++){
                for (int j = 2; j < ny + 2; j++){
                    for (int i = 2; i < nx + 2; i++){
                        const FLOAT dx_0 = _parameters.meshsize->getDx(i  ,j  ,k  );
                        const FLOAT dx_M1= _parameters.meshsize->getDx(i-1,j  ,k  );
                        const FLOAT dx_P1= _parameters.meshsize->getDx(i+1,j  ,k  );
                        const FLOAT dy_0 = _parameters.meshsize->getDy(i  ,j  ,k  );
                        const FLOAT dy_M1= _parameters.meshsize->getDy(i  ,j-1,k  );
                        const FLOAT dy_P1= _parameters.meshsize->getDy(i  ,j+1,k  );
                        const FLOAT dz_0 = _parameters.meshsize->getDz(i  ,j  ,k  );
                        const FLOAT dz_M1= _parameters.meshsize->getDz(i  ,j  ,k-1);
                        const FLOAT dz_P1= _parameters.meshsize->getDz(i  ,j  ,k+1);

                        const FLOAT dx_W = 0.5*(dx_0+dx_M1);
                        const FLOAT dx_E = 0.5*(dx_0+dx_P1);
                        const FLOAT dx_S = 0.5*(dy_0+dy_M1);
                        const FLOAT dx_N = 0.5*(dy_0+dy_P1);
                        const FLOAT dx_B = 0.5*(dz_0+dz_M1);
                        const FLOAT dx_T = 0.5*(dz_0+dz_P1);

                        const FLOAT a_W  =  2.0/(dx_W*(dx_W+dx_E));
                        const FLOAT a_E  =  2.0/(dx_E*(dx_W+dx_E));
                        const FLOAT a_N  =  2.0/(dx_N*(dx_N+dx_S));
                        const FLOAT a_S  =  2.0/(dx_S*(dx_N+dx_S));
                        const FLOAT a_T  =  2.0/(dx_T*(dx_T+dx_B));
                        const FLOAT a_B  =  2.0/(dx_B*(dx_T+dx_B));
                        const FLOAT a_C  = -2.0/(dx_E*dx_W)-2.0/(dx_N*dx_S)-2.0/(dx_B*dx_T);

                        resnorm += pow((
                                     _flowField.getRHS().getScalar(i,j,k)
                                     - a_W*P.getScalar(i-1,j,k) - a_E*P.getScalar(i+1,j,k)
                                     - a_S*P.getScalar(i,j-1,k) - a_N*P.getScalar(i,j+1,k)
                                     - a_B*P.getScalar(i,j,k-1) - a_T*P.getScalar(i,j,k+1)
                                     - a_C*P.getScalar(i,j,k)
                                   ),2);
                    }
                }
            }
            resnorm = sqrt (resnorm / (nx * ny * nz));
            //std::cout << "Residual norm: " << resnorm << std::endl;

            it++;
            iterations--;
        } while ( resnorm > tol && iterations);

    } if (_parameters.geometry.dim == 2){
        do {
            for (int j = 2; j < ny + 2; j++){
                for (int i = 2; i < nx + 2; i++){
                        const FLOAT dx_0 = _parameters.meshsize->getDx(i  ,j  );
                        const FLOAT dx_M1= _parameters.meshsize->getDx(i-1,j  );
                        const FLOAT dx_P1= _parameters.meshsize->getDx(i+1,j  );
                        const FLOAT dy_0 = _parameters.meshsize->getDy(i  ,j  );
                        const FLOAT dy_M1= _parameters.meshsize->getDy(i  ,j-1);
                        const FLOAT dy_P1= _parameters.meshsize->getDy(i  ,j+1);

                        const FLOAT dx_W = 0.5*(dx_0+dx_M1);
                        const FLOAT dx_E = 0.5*(dx_0+dx_P1);
                        const FLOAT dx_S = 0.5*(dy_0+dy_M1);
                        const FLOAT dx_N = 0.5*(dy_0+dy_P1);

                        const FLOAT a_W  =  2.0/(dx_W*(dx_W+dx_E));
                        const FLOAT a_E  =  2.0/(dx_E*(dx_W+dx_E));
                        const FLOAT a_N  =  2.0/(dx_N*(dx_N+dx_S));
                        const FLOAT a_S  =  2.0/(dx_S*(dx_N+dx_S));
                        const FLOAT a_C  = -2.0/(dx_E*dx_W)-2.0/(dx_N*dx_S);

                        const FLOAT gaussSeidel = 1.0/a_C*(_flowField.getRHS().getScalar(i,j)
                                                 - a_W*P.getScalar(i-1,j) - a_E*P.getScalar(i+1,j)
                                                 - a_S*P.getScalar(i,j-1) - a_N*P.getScalar(i,j+1));
                        P.getScalar(i,j) = omg*gaussSeidel + (1.0-omg) * P.getScalar(i,j);
                }
            }

            resnorm = 0.0;
            for (int j = 2; j < ny + 2; j++){
                for (int i = 2; i < nx + 2; i++){
                        const FLOAT dx_0 = _parameters.meshsize->getDx(i  ,j  );
                        const FLOAT dx_M1= _parameters.meshsize->getDx(i-1,j  );
                        const FLOAT dx_P1= _parameters.meshsize->getDx(i+1,j  );
                        const FLOAT dy_0 = _parameters.meshsize->getDy(i  ,j  );
                        const FLOAT dy_M1= _parameters.meshsize->getDy(i  ,j-1);
                        const FLOAT dy_P1= _parameters.meshsize->getDy(i  ,j+1);

                        const FLOAT dx_W = 0.5*(dx_0+dx_M1);
                        const FLOAT dx_E = 0.5*(dx_0+dx_P1);
                        const FLOAT dx_S = 0.5*(dy_0+dy_M1);
                        const FLOAT dx_N = 0.5*(dy_0+dy_P1);

                        const FLOAT a_W  =  2.0/(dx_W*(dx_W+dx_E));
                        const FLOAT a_E  =  2.0/(dx_E*(dx_W+dx_E));
                        const FLOAT a_N  =  2.0/(dx_N*(dx_N+dx_S));
                        const FLOAT a_S  =  2.0/(dx_S*(dx_N+dx_S));
                        const FLOAT a_C  = -2.0/(dx_E*dx_W)-2.0/(dx_N*dx_S);

                        const FLOAT residual = _flowField.getRHS().getScalar(i,j)
                                               - a_W*P.getScalar(i-1,j) - a_E*P.getScalar(i+1,j)
                                               - a_S*P.getScalar(i,j-1) - a_N*P.getScalar(i,j+1)
                                               - a_C*P.getScalar(i,j);
                        resnorm += residual*residual;
                }
            }
            resnorm = sqrt (resnorm / (nx * ny));

            for ( int j = 2; j < ny + 2; j++ ) {
                P.getScalar(1,j) = P.getScalar(2,j);
                P.getScalar(nx+2,j) = P.getScalar(nx+1,j);
            }

            for ( int i = 2; i < nx + 2; i++ ) {
                P.getScalar(i,1) = P.getScalar(i,2);
                P.getScalar(i,ny+2) = P.getScalar(i,ny+1);
            }

            //std::cout << "Residual norm: " << resnorm << std::endl;
            iterations--;
            it++;
        } while ( resnorm > tol && iterations);
    }
}
