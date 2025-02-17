// MIT License
// 
// Copyright (c) 2025 molecularinformatics  
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


#include <stdio.h>
#include "cuda_functions.cuh"



////////////////////////////////////////////////////////////////////////////////
/// Math helper functions
////////////////////////////////////////////////////////////////////////////////


__device__ void matvec3x3x3(float mat[3][3], const float * vec, float * result){

    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = 0.0;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i] += mat[i][j] * vec[j];
        }
    }
}

__device__ void cross_product(const float a[3], const float b[3], float c[3]){
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}


__device__ void quaternion_to_rotation_matrix(const float * q, float M[3][3]){

    // temp variables to make more readable
    float w = q[0];
    float x = q[1];
    float y = q[2];
    float z = q[3];

    // Compute rotation matrix elements
    M[0][0] = 1 - 2*y*y - 2*z*z;
    M[0][1] = 2*x*y - 2*w*z;
    M[0][2] = 2*x*z + 2*w*y;
    M[1][0] = 2*x*y + 2*w*z;
    M[1][1] = 1 - 2*x*x - 2*z*z;
    M[1][2] = 2*y*z - 2*w*x;
    M[2][0] = 2*x*z - 2*w*y;
    M[2][1] = 2*y*z + 2*w*x;
    M[2][2] = 1 - 2*x*x - 2*y*y;
}

__device__ void transform_inplace(float mat[3][3], float * mol, int Nmol){
    // only transforms the first 3 values in the array! The 4th value is not a 
    // coordinate so we do not transform it
    const int D=4;

    for(int i=0;i<Nmol;i++){

        float temp[3] = {mol[i*D], mol[i*D+1], mol[i*D+2]};
        matvec3x3x3(mat, &mol[i*D], temp);
        mol[i*D]   = temp[0];
        mol[i*D+1] = temp[1];
        mol[i*D+2] = temp[2];
    }
}

__device__ void translate_inplace(float t[3], float * mol, int Nmol){
    // only transforms the first 3 values in the array! The 4th value is not a 
    // coordinate so we do not transform it
    const int D=4;

    for(int i=0;i<Nmol;i++){

        mol[i*D]   += t[0];
        mol[i*D+1] += t[1];
        mol[i*D+2] += t[2];

    }
}

__device__ void axis_angle_to_quat(const float * axis_angle, float * q){

    // axis must be normalized. We assume it is
    
    float hangle = axis_angle[3]*0.5;
    q[0] = cos(hangle);
    q[1] = axis_angle[0]*sin(hangle);
    q[2] = axis_angle[1]*sin(hangle);
    q[3] = axis_angle[2]*sin(hangle);

}




__device__ void start_mode_transform(const float * mol, float * transformed_mol, float * transform, int Nmol, int method, int index){

    // set transform to default values
    transform[0]=1;
    transform[1]=0;
    transform[2]=0;
    transform[3]=0;
    transform[4]=0;
    transform[5]=0;
    transform[6]=0;


    // initialise the transformed one to be the same as the original one
    for(int i=0; i<Nmol; i++){
        transformed_mol[D*i] = mol[D*i];
        transformed_mol[D*i+1] = mol[D*i+1];
        transformed_mol[D*i+2] = mol[D*i+2];
        transformed_mol[D*i+3] = mol[D*i+3];
    }

    // number of start mode configuration for the chosen method
    //int method_n = start_mode_n[method];

    // the methods have a different number of start positions
    // they are listed in here as consts
    // TODO: can hardcode the qs and matrices
    switch (method){

        case 0:
        {
            // method 0
            // do nothing
            break;
        }
        
        case 1:
        {
            // list the transformations
            const float transforms_1[16] = {
                1,0,0, 0,
                1,0,0,PI,
                0,1,0,PI,
                0,0,1,PI
            };

            // get the one we need
            const float * my_transform = &transforms_1[index*4];

            // apply the transformation
            float q[4];
            axis_angle_to_quat(my_transform,q);

            float M[3][3];
            quaternion_to_rotation_matrix(q,M);

            for(int i=0; i<Nmol; i++){
                matvec3x3x3(M, &mol[D*i], &transformed_mol[D*i]);
            }

            transform[0] = q[0];
            transform[1] = q[1];
            transform[2] = q[2];
            transform[3] = q[3];

            break;
        }

        case 2:
        {
  
            // list the transformations
            const float transforms_2[40] = {
                1,0,0,0,
                1,0,0,PI,
                0,1,0,PI,
                0,0,1,PI,
                1,0,0,0.5f*PI,
                0,1,0,0.5f*PI,
                0,0,1,0.5f*PI,
                1,0,0,-0.5f*PI,
                0,1,0,-0.5f*PI,
                0,0,1,-0.5f*PI
            };


            // get the one we need
            const float * my_transform = &transforms_2[index*4];

            // apply the transformation
            float q[4];
            axis_angle_to_quat(my_transform,q);

            float M[3][3];
            quaternion_to_rotation_matrix(q,M);

            for(int i=0; i<Nmol; i++){
                matvec3x3x3(M, &mol[D*i], &transformed_mol[D*i]);
            }

            transform[0] = q[0];
            transform[1] = q[1];
            transform[2] = q[2];
            transform[3] = q[3];

            break;
        }
        default:
            break;

    } // switch
}



///////////////////////////////////////////////////////////////////////////////
/// Volume functions
///////////////////////////////////////////////////////////////////////////////


// shape and color overap of molA and molB
__device__ void volume_single(const float * molA, const int * molA_type, int molA_num_atoms, int NmolA, 
                              const float * molB, const int * molB_type, int molB_num_atoms, int NmolB,
                              const float * rmat, const float * pmat, int N_features, float * v){

    const float PI = 3.14159265358;
    const float KAPPA = 2.41798793102;
    const float CARBONRADII2 = 1.7*1.7;
    const float A = KAPPA/CARBONRADII2;
    const float CONSTANT = pow(PI/(2*A), 1.5);


    float vs = 0.0;
    float vc = 0.0;
    
    // normal atoms first then color
    for(int j=0; j<molB_num_atoms; j++){
        for(int i=0; i<molA_num_atoms; i++){

            float dx = molA[i*D]   - molB[j*D];
            float dy = molA[i*D+1] - molB[j*D+1];
            float dz = molA[i*D+2] - molB[j*D+2];

            float d2 = dx*dx + dy*dy + dz*dz;        

            auto a1 = A;
            auto a2 = A;

            float wa = molA[i*D+3];
            float wb = molB[j*D+3];

            float kij = exp(-a1*a2*d2/(a1+a2))*wa*wb;

            float vij = 8*kij*CONSTANT;

            vs += vij;
        }
    }

    for(int j=molB_num_atoms; j<NmolB; j++){
        for(int i=molA_num_atoms; i<NmolA; i++){

            float dx = molA[i*D]   - molB[j*D];
            float dy = molA[i*D+1] - molB[j*D+1];
            float dz = molA[i*D+2] - molB[j*D+2];

            float d2 = dx*dx + dy*dy + dz*dz;        

            int ta = molA_type[i];
            int tb = molB_type[j];

            float a = rmat[ta*N_features+tb];
            float p = pmat[ta*N_features+tb];
        
            float wa = molA[i*D+3];
            float wb = molB[j*D+3];

            float kij = exp(-a*a*d2/(a+a))*wa*wb;

            float constant = pow(PI/(2*a), 1.5);

            float vij = p*p*kij*constant;

            vc += vij;
        }
    }

    v[0] = vs;
    v[1] = vc;
}




////////////////////////////////////////////////////////////////////////////////
/// gradient functions
////////////////////////////////////////////////////////////////////////////////



// compute shape and color overlap gradients of molA and (molB transformed by q+t) 
// w.r.t q and t
__device__ void get_gradient(const float * molA, const int * molA_type, int molA_num_atoms,int NmolA, 
                             const float * molB, const int * molB_type, int molB_num_atoms,int NmolB,
                             const float * rmat, const float * pmat, int N_features, 
                             float * gradient, float * gradient_color, const float * q, const float * t, float * v){
    

    const float PI = 3.14159265358;
    const float KAPPA = 2.41798793102;
    const float CARBONRADII2 = 1.7*1.7;
    const float A = KAPPA/CARBONRADII2;
    const float CONSTANT = pow(PI/(2*A), 1.5);  

    float vs = 0.0;
    float vc = 0.0;


    float M[3][3] = {{0.0,0.0,0.0},
                     {0.0,0.0,0.0},
                     {0.0,0.0,0.0}};

    quaternion_to_rotation_matrix(q, M);

    #ifdef DEBUG
        printf("loop limits1: %d %d\n",molA_num_atoms, molB_num_atoms);
    #endif

    // loop over molB in outer loop to minimize global memory access
    for( int j=0; j<molB_num_atoms; j++){

        // store values in registers
        float molBTr[3];

        // molB needs to be transformed by q
        matvec3x3x3(M, &molB[j*D], molBTr); 

        float wb = molB[j*D+3];

        // loop over molA in inner loop. Note it is in shared memory
        for(int i=0; i <molA_num_atoms; i++){

           
            // molA needs to be translated by -t
            float molATx = molA[i*D]   - t[0];
            float molATy = molA[i*D+1] - t[1];
            float molATz = molA[i*D+2] - t[2];

            float wa = molA[i*D+3];

            float dx = molATx - molBTr[0];
            float dy = molATy - molBTr[1];
            float dz = molATz - molBTr[2];

            float d2 = dx*dx + dy*dy + dz*dz;
            
        
            // float a1 = A;
            // float a2 = A;
            const float factor = A*A/(A+A);
            
            // float kij = __expf(-a1*a2*d2/(a1+a2))*wa*wb;
            float kij = __expf(-d2*factor)*wa*wb;
            
            float vij = 8*kij*CONSTANT;

            // cross product method to get gradient

            float x = molBTr[0];
            float y = molBTr[1];
            float z = molBTr[2];

            float C[3] = {2*x,2*y,2*z};
            float mv[3] = {0.0,0.0,0.0};
            float delta[3] = {dx, dy, dz};

            cross_product(C, delta, mv);

            // float temp = -2.0*(a1*a2)/(a1+a2)*vij;
            float temp = -2.0*factor*vij;

            // dq[0] is always zero
            // dq[1-3]: 
            gradient[1] += temp*mv[0];
            gradient[2] += temp*mv[1];
            gradient[3] += temp*mv[2];

            // d[x,y,z]:
            gradient[4] += temp*delta[0];
            gradient[5] += temp*delta[1];
            gradient[6] += temp*delta[2];

            vs += vij;
        }
    }

    #ifdef DEBUG
    printf("loop limits2: %d:%d %d:%d\n",molA_num_atoms, NmolA, molB_num_atoms, NmolB);
    #endif
    // TODO: combine both loops to half the global memory access
    for( int j=molB_num_atoms; j<NmolB; j++){

        // store values in registers
        float molBTr[3];

        // molB needs to be transformed by q
        matvec3x3x3(M, &molB[j*D], molBTr); 

        float wb = molB[j*D+3];

        // TODO: can terminate j loop early
        int tb = molB_type[j];


        // loop over molA in inner loop. Note it is in shared memory
        for(int i=molA_num_atoms; i<NmolA; i++){

            int ta = molA_type[i];

            float a = rmat[ta*N_features+tb];
            float p = pmat[ta*N_features+tb];

            // if (a==0.0f) continue; //TODO


            // molA needs to be translated by -t
            float molATx = molA[i*D]   - t[0];
            float molATy = molA[i*D+1] - t[1];
            float molATz = molA[i*D+2] - t[2];

            float wa = molA[i*D+3];


            float dx = molATx - molBTr[0];
            float dy = molATy - molBTr[1];
            float dz = molATz - molBTr[2];

            float d2 = dx*dx + dy*dy + dz*dz;
            



            float inv2a = 1.0/(2*a);
            float factor = a*a*inv2a;
        
            float kij = __expf(-factor*d2)*wa*wb;

            // float constant = pow(PI/(2*a), 1.5);

            // this should be faster
            // float constant = __powf(PI/(2*a), 1.5);

            // this should be even faster
            float temp1 = PI*inv2a;
            float constant = temp1*__fsqrt_rn(temp1);
            // TODO: should just precompute all factors and constants for all rmat 

            float vij = p*p*kij*constant;


            float x = molBTr[0];
            float y = molBTr[1];
            float z = molBTr[2];

            float C[3] = {2*x,2*y,2*z};
            float mv[3] = {0.0,0.0,0.0};
            float delta[3] = {dx, dy, dz};

            cross_product(C, delta, mv);

            // float temp = -2.0*(a*a)/(a+a)*vij;
            float temp = -2.0*factor*vij;

            // dq[0] is always zero
            // dq[1-3]: 
            gradient_color[1] += temp*mv[0];
            gradient_color[2] += temp*mv[1];
            gradient_color[3] += temp*mv[2];

            // d[x,y,z]:
            gradient_color[4] += temp*delta[0];
            gradient_color[5] += temp*delta[1];
            gradient_color[6] += temp*delta[2];

            vc += vij;
        }
    }

    v[0] = vs;
    v[1] = vc; 

}


////////////////////////////////////////////////////////////////////////////////
/// Optimization functions
////////////////////////////////////////////////////////////////////////////////

// Adagrad optimization step
__device__ void adagrad_step(float *q, float *t,  float * g, 
                  float * cache, float lr_q, float lr_t) {

    
    for(int i=0; i<7; ++i){
        cache[i] =  cache[i] + g[i]*g[i];
    }

    q[0] -= lr_q*g[0]/(sqrt(cache[0])+EPSILON);
    q[1] -= lr_q*g[1]/(sqrt(cache[1])+EPSILON);
    q[2] -= lr_q*g[2]/(sqrt(cache[2])+EPSILON);
    q[3] -= lr_q*g[3]/(sqrt(cache[3])+EPSILON);

    t[0] -= lr_t*g[4]/(sqrt(cache[4])+EPSILON);
    t[1] -= lr_t*g[5]/(sqrt(cache[5])+EPSILON);
    t[2] -= lr_t*g[6]/(sqrt(cache[6])+EPSILON);
}

// optimization kernel. Entry point to the cude kernel code from the c++  cuda
// interface code. 
__global__ void optimize(const float * molA_global, const int * molA_type_global, const int molA_num_atoms_i,  int NmolA, 
                         const float * molBs,       const int * molB_types,       const int * molB_num_atoms,  int NmolB, 
                         const float * rmat_global, const float * pmat_global, int N_features,
                         float * scores, int Nv, bool optim_color, float lr_q, float lr_t, int nsteps, float mixing_param, int start_mode){
    

    // each thead works on one overlap
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

  
    // TODO compile time constant
    const int DV = 20;

    // shared memory for MolA
    extern __shared__ float shared[];

    // can only have one shared so we need to store multiple arrays in it
    float * molA = shared;
    float * qr = (float *)&shared[NmolA*D];    
    float * rmat = (float *)&shared[NmolA*D+8];
    float * pmat = (float *)&shared[NmolA*D+8+N_features*N_features];
    int * molA_type = (int *)&shared[NmolA*D+8+2*N_features*N_features];

    size_t tidx = threadIdx.x;

    // printf("blockIdx.x: %d, blockDim.x: %d, threadIdx.x: %d\n", 
        //    blockIdx.x, blockDim.x, threadIdx.x);

    if (tidx < NmolA){
        // copy molA to shared memory
        //for(int l=tidx; l<NmolA*D; l+=blockDim.x){
        //    molA[l] = molA_global[l];
        //} // done by the start mode part now
        
        // copy molA_type to shared memory
        for(int l=tidx; l<NmolA; l+=blockDim.x){
            molA_type[l] = molA_type_global[l];
        }
    }

    if (tidx < (N_features*N_features)){
        // copy rmat and pmat to shared memory
        for(int l=tidx; l<(N_features*N_features); l+=blockDim.x){
            rmat[l] = rmat_global[l];
            pmat[l] = pmat_global[l];
        }
    }

    __syncthreads();
    

    // loop over the start modes

    int N_starts = start_mode_n[start_mode];

    // register for scores
    float scores_register[20];
    for(int j=0;j<20;++j){
        scores_register[j]=0.0f;
    }
    

    for(int n=0; n<N_starts; n++){
        // create the new A and record the transformation
        // only one thread needs to do this #TODO can be sped up by splitting over threads
        if (tidx==0){
            start_mode_transform(molA_global, molA, qr, NmolA, start_mode, n);
        }
        __syncthreads();

        #ifdef DEBUG
        printf("in __optimize__, start mode index %d / %d\n", n, N_starts);
        #endif
        if (idx < Nv){

            const float * molB = &molBs[idx*NmolB*D];
            const int * molB_type = &molB_types[idx*NmolB];
            int molB_num_atom = molB_num_atoms[idx];


            // self overlaps
            float va[2];
            float vb[2];

            volume_single(molA, molA_type, molA_num_atoms_i, NmolA,
                        molA, molA_type, molA_num_atoms_i, NmolA,
                        rmat, pmat, N_features, va);

            volume_single(molB, molB_type, molB_num_atom, NmolB,
                        molB, molB_type, molB_num_atom, NmolB,
                        rmat, pmat, N_features, vb);
            

            #ifdef DEBUG
            printf("idx: %d\n", idx);
            #endif

            float q[4] = {1.0,0.0,0.0,0.0}; // initial q
            float t[3] = {0.0,0.0,0.0}; // initial t




            #ifdef DEBUG
            printf("molB_num_atom %d\n", molB_num_atom);
            #endif

            // pair overlap volumes [shape, color]
            float vs[2] = {0.0,0.0};

            float cache[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0}; //for optimizer

            float g_factor   = 1.0/(va[0]+vb[0]);

            // this can NaN if the self color is zero
            float g_c_factor;
            if(optim_color) 
                g_c_factor = 1.0/(va[1]+vb[1]);
            else
                g_c_factor = 0.0;

            for(int k=0; k < nsteps; ++k){

                float g[7]   = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
                float g_c[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
                float g_combo[7] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};

                get_gradient(molA, molA_type, molA_num_atoms_i, NmolA, 
                            molB, molB_type, molB_num_atom,    NmolB, 
                            rmat, pmat, N_features,
                            g, g_c, q, t, vs);

                // 3. update
                
                // normalize and combine

                #pragma unroll
                for(int i=0;i<7;++i){
                    g[i] = g[i]*g_factor;
                }
                
                if (optim_color){
                
                    #pragma unroll
                    for(int i=0;i<7;++i){
                        g_c[i] = g_c[i]*g_c_factor;
                    }
                    
                    
                    #pragma unroll
                    for(int i=0;i<7;++i){
                        g_combo[i] = (1-mixing_param)*g[i]+ mixing_param*g_c[i];
                    }
                    
                    
                    adagrad_step(q,t,g_combo, cache, lr_q, lr_t);
                }
                else{
                    adagrad_step(q,t,g, cache, lr_q, lr_t);
                }

                
                //normalize q
                float magq = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] * q[3]*q[3]);

                // std::cout << magq << std::endl;
                q[0] = q[0]/magq;
                q[1] = q[1]/magq;
                q[2] = q[2]/magq;
                q[3] = q[3]/magq;

            }
            // final step to get latest volumes
            // get_gradient computes the volumes, we dont care about the gradients it returns
            float temp[7]   = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
            get_gradient(molA, molA_type, molA_num_atoms_i, NmolA, 
                            molB, molB_type, molB_num_atom, NmolB, 
                            rmat, pmat, N_features,
                            temp,temp, q, t, vs);


            
            // compute tamimoto scores

            // shape
            float ts = vs[0] / (va[0]+vb[0] - vs[0]);

            // color
            float tc = 0.0;

            if(optim_color){
                tc = vs[1] / (va[1] + vb[1] - vs[1]);
            }

            // we use the mixing param to weight the tanimotos
            // this is the objective function we have optimized
            float total = (ts *(1-mixing_param) + tc*mixing_param);


            #ifdef DEBUG
            printf("completed optim on %d\n", idx);
            #endif

            // check if this start mode is the best
            if (total > scores_register[0]){
                // store the scores, volumes, and transformation
                scores_register[0] = total; // combination tanimoto of shape and color
                scores_register[1] = ts; // shape tanimoto
                scores_register[2] = tc; // color tanimoto
                scores_register[3] = vs[0]; // volume shape
                scores_register[4] = vs[1]; // volumes color
                scores_register[5] = va[0]; // self a
                scores_register[6] = vb[0]; // self b
                scores_register[7] = va[1]; // self a color
                scores_register[8] = vb[1]; // self b color
                scores_register[9]  = q[0]; // fit q
                scores_register[10] = q[1];
                scores_register[11] = q[2];
                scores_register[12] = q[3];
                scores_register[13] = t[0]; // fit t
                scores_register[14] = t[1];
                scores_register[15] = t[2];
                scores_register[16] = qr[0]; // start_mode q
                scores_register[17] = qr[1];
                scores_register[18] = qr[2];
                scores_register[19] = qr[3];



                #ifdef DEBUG
                printf("volumes on %d: %f %f %f %f %f %f\n", idx, vs[0], vs[1], va[0], vb[0], va[1], vb[1]);
                #endif
            }

        

        scores[idx*DV]    = scores_register[0];   
        scores[idx*DV+1]  = scores_register[1]; 
        scores[idx*DV+2]  = scores_register[2]; 
        scores[idx*DV+3]  = scores_register[3]; 
        scores[idx*DV+4]  = scores_register[4]; 
        scores[idx*DV+5]  = scores_register[5]; 
        scores[idx*DV+6]  = scores_register[6]; 
        scores[idx*DV+7]  = scores_register[7]; 
        scores[idx*DV+8]  = scores_register[8]; 
        scores[idx*DV+9]  = scores_register[9]; 
        scores[idx*DV+10] = scores_register[10];
        scores[idx*DV+11] = scores_register[11];
        scores[idx*DV+12] = scores_register[12];
        scores[idx*DV+13] = scores_register[13];
        scores[idx*DV+14] = scores_register[14];
        scores[idx*DV+15] = scores_register[15];
        scores[idx*DV+16] = scores_register[16];
        scores[idx*DV+17] = scores_register[17];
        scores[idx*DV+18] = scores_register[18];
        scores[idx*DV+19] = scores_register[19];

        }

    __syncthreads();
    }
}






// This is the entry point function. The optimization kernel is lauched from here
void optimize_overlap_gpu(const float * molA, const int * molA_type, int molA_num_atoms,  int NmolA, 
                          const float * molBs,const int * molB_types, const int * molB_num_atoms, int NmolB, long num_molBs, 
                          const float * rmat, const float * pmat, int N_features, float * scores, bool optim_color, float lr_q, float lr_t, int nsteps, float mixing_param, int start_mode, int device_id){


    #ifdef DEBUG
    std::cout << " in optimize_overlap_gpu " << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    #endif


    long Nv = num_molBs;

    dim3 blockDim(NTHREADS); // threads per block
    dim3 gridDim((Nv + blockDim.x-1)/blockDim.x); // number of blocks

    const int D = 4; // TODO: compile time constant
    size_t sharedMemSize = NmolA*D*sizeof(float)+8*sizeof(float)+NmolA*sizeof(int)+2*(N_features*N_features)*sizeof(float);
   
    //TODO: should check this is not too large
    // std::cout << " sharedMemSize = " << sharedMemSize << " at " << __FILE__ << ":" << __LINE__ << std::endl;
    
    cudaSetDevice(device_id); optimize<<<gridDim,blockDim,sharedMemSize>>>(molA,  molA_type,  molA_num_atoms, NmolA,
                                                 molBs, molB_types, molB_num_atoms,  NmolB, 
                                                 rmat,  pmat, N_features,
                                                 scores, Nv, optim_color,lr_q, lr_t, nsteps, mixing_param, start_mode);
    cudaSetDevice(device_id); cudaGetLastError();
    cudaSetDevice(device_id); cudaDeviceSynchronize();
    cudaSetDevice(device_id); cudaGetLastError();

}


// // wrapper function for testing
// void get_gradient_gpu(const float * molA, const int * molA_type, int molA_num_atoms,  int NmolA, 
//                       const float * molB, const int * molB_type, int molB_num_atoms,  int NmolB, 
//                       const float * q, const float * t, float * gradient, float * vs){
    

//     get_gradient_test_kernel_wrapper<<<1,1>>>(molA, molA_type, molA_num_atoms, NmolA, 
//                                        molB, molB_type, molB_num_atoms, NmolB, 
//                                        q, t, gradient, vs);
//     CUDA_CHECK_ERROR(cudaGetLastError());
//     cudaDeviceSynchronize();
//     CUDA_CHECK_ERROR(cudaGetLastError());

// }
