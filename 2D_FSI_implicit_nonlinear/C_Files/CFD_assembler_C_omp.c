/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch>
 */

#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>
#include "Tools.h"

#define INVJAC(i,j,k) invjac[i+(j+k*dim)*noe]
#define GRADREFPHIV(i,j,k) gradrefphiV[i+(j+k*NumQuadPoints)*nlnV]
#define GRADREFPHIP(i,j,k) gradrefphiP[i+(j+k*NumQuadPoints)*nlnP]
#ifdef _OPENMP
#include <omp.h>
#else
#warning "OpenMP not enabled. Compile with mex CFD_assembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp""
#endif


/*************************************************************************/
void AssembleStokes(mxArray* plhs[], const mxArray* prhs[])
{
    double* dim_ptr = mxGetPr(prhs[2]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[3]);
    double* nln_ptrV = mxGetPr(prhs[4]);
    int nlnV     = (int)(nln_ptrV[0]);
    double* nln_ptrP = mxGetPr(prhs[5]);
    int nlnP     = (int)(nln_ptrP[0]);
    int numRowsElements  = mxGetM(prhs[3]);
    
    int nln  = nlnV + nlnP;
    int nln2 = nln*nln;
    int local_matrix_size = nlnV*nlnV*dim*dim + 2*nlnV*nlnP*dim;
    int global_lenght = noe * local_matrix_size;
    
    plhs[0] = mxCreateDoubleMatrix(global_lenght,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(global_lenght,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(global_lenght,1, mxREAL);
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);
    
    int NumQuadPoints     = mxGetN(prhs[7]);
    
    double* NumNodes_ptr = mxGetPr(prhs[6]);
    int NumScalarDofsV     = (int)(NumNodes_ptr[0] / dim);
    
    double* w   = mxGetPr(prhs[7]);
    double* invjac = mxGetPr(prhs[8]);
    double* detjac = mxGetPr(prhs[9]);
    double* phiV = mxGetPr(prhs[10]);
    double* gradrefphiV = mxGetPr(prhs[11]);
    double* phiP = mxGetPr(prhs[12]);
    
    double* elements  = mxGetPr(prhs[3]);
        
    double* material_param = mxGetPr(prhs[1]);
    double viscosity = material_param[0];
    
    /* Assembly: loop over the elements */
    int ie;
    
#pragma omp parallel for shared(invjac,detjac,elements,myAcols,myArows,myAcoef) private(ie) firstprivate(phiV,phiP,gradrefphiV,w,numRowsElements,local_matrix_size,nlnV,nlnP,NumQuadPoints,NumScalarDofsV,viscosity,dim)
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        int k, q, d1, d2;
        
        double gradphiV[NumQuadPoints][nlnV][dim];
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            for (k = 0; k < nlnV; k = k + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphiV[q][k][d1] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphiV[q][k][d1] = gradphiV[q][k][d1] + INVJAC(ie,d1,d2)*GRADREFPHIV(k,q,d2);
                    }
                }
            }
        }
        
        int iii = 0;
        int ii = 0;
        int a, b;
        
        /* loop over velocity test functions --> a */
        for (a = 0; a < nlnV; a = a + 1 )
        {
            /* loop over velocity trial functions --> b */
            for (b = 0; b < nlnV; b = b + 1 )
            {
                double aloc[dim][dim];
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        aloc[d1][d2] = 0;
                    }
                }
                for (q = 0; q < NumQuadPoints; q = q + 1 )
                {
                    double scalar_lapl = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        scalar_lapl  = scalar_lapl + viscosity * gradphiV[q][a][d2] * gradphiV[q][b][d2]  * w[q];
                    }
                    
                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                        {
                            aloc[d1][d2]  = aloc[d1][d2] + viscosity * gradphiV[q][a][d2] * gradphiV[q][b][d1]  * w[q];
                        }
                    }
                    
                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        aloc[d1][d1]  = aloc[d1][d1] + scalar_lapl;
                    }
                    
                }
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        myArows[ie*local_matrix_size+iii] = elements[a+ie*numRowsElements] + d1 * NumScalarDofsV;
                        myAcols[ie*local_matrix_size+iii] = elements[b+ie*numRowsElements] + d2 * NumScalarDofsV;
                        myAcoef[ie*local_matrix_size+iii] = aloc[d1][d2]*detjac[ie];
                        
                        iii = iii + 1;
                    }
                }
            }
            
            /* loop over pressure trial functions --> b */
            for (b = 0; b < nlnP; b = b + 1 )
            {
                double alocDiv[dim];
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    alocDiv[d1] = 0;
                }
                
                for (q = 0; q < NumQuadPoints; q = q + 1 )
                {
                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        alocDiv[d1]  = alocDiv[d1] + ( - phiP[b+q*nlnP] * gradphiV[q][a][d1] ) * w[q];
                    }
                }
                
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    myArows[ie*local_matrix_size+iii] = elements[a+ie*numRowsElements] + d1 * NumScalarDofsV;
                    myAcols[ie*local_matrix_size+iii] = elements[b+ie*numRowsElements] + dim * NumScalarDofsV;
                    myAcoef[ie*local_matrix_size+iii] = alocDiv[d1]*detjac[ie];
                    
                    iii = iii + 1;
                    
                    myAcols[ie*local_matrix_size+iii] = elements[a+ie*numRowsElements] + d1  * NumScalarDofsV;
                    myArows[ie*local_matrix_size+iii] = elements[b+ie*numRowsElements] + dim * NumScalarDofsV;
                    myAcoef[ie*local_matrix_size+iii] = -alocDiv[d1]*detjac[ie];
                    
                    iii = iii + 1;
                }
                
            }
        }
    }
}

/*************************************************************************/
void AssembleConvectiveALE(mxArray* plhs[], const mxArray* prhs[])
{
    double* dim_ptr = mxGetPr(prhs[2]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[3]);
    double* nln_ptrV = mxGetPr(prhs[4]);
    int nlnV     = (int)(nln_ptrV[0]);
    int numRowsElements  = mxGetM(prhs[3]);
    
    int local_matrix_size1 = nlnV*nlnV*dim;
    int local_matrix_size2 = nlnV*nlnV*dim*dim;

    int global_lenght1 = noe * local_matrix_size1;
    int global_lenght2 = noe * local_matrix_size2;

    plhs[0] = mxCreateDoubleMatrix(global_lenght1,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(global_lenght1,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(global_lenght1,1, mxREAL);
    
    plhs[3] = mxCreateDoubleMatrix(global_lenght2,1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(global_lenght2,1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(global_lenght2,1, mxREAL);
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);
    
    double* myBrows    = mxGetPr(plhs[3]);
    double* myBcols    = mxGetPr(plhs[4]);
    double* myBcoef    = mxGetPr(plhs[5]);
    
    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    
    double* NumNodes_ptr = mxGetPr(prhs[5]);
    int NumScalarDofsV     = (int)(NumNodes_ptr[0] / dim);
    
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phiV = mxGetPr(prhs[9]);
    double* gradrefphiV = mxGetPr(prhs[10]);
    double* U_h   = mxGetPr(prhs[11]);
    double* Conv_velocity   = mxGetPr(prhs[12]);
    
    double* elements  = mxGetPr(prhs[3]);
    
    double* material_param = mxGetPr(prhs[1]);
    double density = material_param[0];
    
    /* Assembly: loop over the elements */
    int ie, d1, d2;
        
#pragma omp parallel for shared(invjac,detjac,elements,myAcols,myArows,myAcoef,myBrows,myBcols,myBcoef,U_h,Conv_velocity) private(ie,k,l,q,d1,d2) firstprivate(phiV,gradrefphiV,w,numRowsElements,local_matrix_size1,local_matrix_size2,nlnV,NumScalarDofsV,density)
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        double gradphiV[dim][nlnV][NumQuadPoints];
        double U_hq[dim][NumQuadPoints];
        double GradV[dim][dim];
        double GradU[dim][dim];
        double GradUh[dim][dim][NumQuadPoints];
        double ConvVel_hq[dim][NumQuadPoints];

        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* Compute Gradient of Basis functions*/
            for (k = 0; k < nlnV; k = k + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphiV[d1][k][q] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphiV[d1][k][q] = gradphiV[d1][k][q] + INVJAC(ie,d1,d2)*GRADREFPHIV(k,q,d2);
                    }
                }
            }
            /* Compute U_h and Grad(U_h) on the quadrature nodes of the current element*/
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                U_hq[d1][q] = 0;
                ConvVel_hq[d1][q] = 0;
                for (k = 0; k < nlnV; k = k + 1 )
                {
                    int e_k;
                    e_k = (int)(elements[ie*numRowsElements + k] + d1*NumScalarDofsV - 1);
                    U_hq[d1][q] = U_hq[d1][q] + U_h[e_k] * phiV[k+q*nlnV];
                    ConvVel_hq[d1][q] = ConvVel_hq[d1][q] + Conv_velocity[e_k] * phiV[k+q*nlnV];
                }
                
                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                {
                    GradUh[d1][d2][q] = 0;
                    for (k = 0; k < nlnV; k = k + 1 )
                    {
                        int e_k;
                        e_k = (int)(elements[ie*numRowsElements + k] + d1*NumScalarDofsV - 1);
                        GradUh[d1][d2][q] = GradUh[d1][d2][q] + U_h[e_k] * gradphiV[d2][k][q];
                    }
                }
            }
        }
        
        int iii = 0, iii2 = 0;
        int a, b, i_c, j_c;
        
        /* loop over velocity test functions --> a */
        for (a = 0; a < nlnV; a = a + 1 )
        {
            /* Assemble C1 = U_h * grad(u) * v */
            /* loop over velocity trial functions --> b */
            for (b = 0; b < nlnV; b = b + 1 )
            {
                double aloc = 0;
                for (q = 0; q < NumQuadPoints; q = q + 1 )
                {
                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        aloc  = aloc + ConvVel_hq[d1][q] * gradphiV[d1][b][q] * phiV[a+q*nlnV] * w[q];
                    }
                }
                
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    myArows[ie*local_matrix_size1+iii] = elements[a+ie*numRowsElements] + d1 * NumScalarDofsV;
                    myAcols[ie*local_matrix_size1+iii] = elements[b+ie*numRowsElements] + d1 * NumScalarDofsV;
                    myAcoef[ie*local_matrix_size1+iii] = density*aloc*detjac[ie];
                    
                    iii = iii + 1;
                }
            }
            
            /* Assemble C2 = u * grad(U_h) * v */
            /* loop over test velocity components --> i_c */
            for (i_c = 0; i_c < dim; i_c = i_c + 1 )
            {
                /* loop over velocity trial functions --> b */
                for (b = 0; b < nlnV; b = b + 1 )
                {
                    /* loop over test velocity components --> j_c */
                    for (j_c = 0; j_c < dim; j_c = j_c + 1 )
                    {
                        double aloc = 0;
                        for (q = 0; q < NumQuadPoints; q = q + 1 )
                        {
                            aloc  = aloc + phiV[b+q*nlnV] * GradUh[i_c][j_c][q] * phiV[a+q*nlnV] * w[q];
                        }
                        myBrows[ie*local_matrix_size2+iii2] = elements[a+ie*numRowsElements] + i_c * NumScalarDofsV;
                        myBcols[ie*local_matrix_size2+iii2] = elements[b+ie*numRowsElements] + j_c * NumScalarDofsV;
                        myBcoef[ie*local_matrix_size2+iii2] = density*aloc*detjac[ie];
                        
                        iii2 = iii2 + 1;
                        
                    }
                }
            }
                
        }
    }
}

/*************************************************************************/
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    char *Assembly_name = mxArrayToString(prhs[0]);
            
    if (strcmp(Assembly_name, "Stokes")==0)
    {
        /* Check for proper number of arguments */
        if(nrhs!=13) {
            mexErrMsgTxt("13 inputs are required.");
        } else if(nlhs>3) {
            mexErrMsgTxt("Too many output arguments.");
        }

        AssembleStokes(plhs, prhs);
    }  
            
    if (strcmp(Assembly_name, "convectiveALE")==0)
    {
        /* Check for proper number of arguments */
        if(nrhs!=13) {
            mexErrMsgTxt("13 inputs are required.");
        } else if(nlhs>6) {
            mexErrMsgTxt("Too many output arguments.");
        }
        
        AssembleConvectiveALE(plhs, prhs);
    }  

    
    mxFree(Assembly_name);
}
/*************************************************************************/

