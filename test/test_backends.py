import numpy as np
import torch
from _roshambo2_cpp  import test_overlap as overlap_cpp, test_gradient as test_gradient_cpu
#from _roshambo2_cuda import overlap as overlap_cuda, test_gradient as test_gradient_gpu
from roshambo2.backends._pytorch_backend import overlap_volume, quaternion_to_rotation_matrix
import time
import copy

NPDTYPE_CPP = np.float32
NPDTYPE_CUDA = np.float32
TORCHDTYPE = torch.float32


# test volume between cpp and pytorch
def test_vols():
    Nmols=10000
    NB = 30
    a = np.random.rand(1,20,3)*10
    b = np.random.rand(Nmols,NB,3)*10
    w = np.ones((Nmols,NB))

    a = a.astype(dtype=NPDTYPE_CPP)
    b = b.astype(dtype=NPDTYPE_CPP)
    w = w.astype(dtype=NPDTYPE_CPP)

    print(b.shape)

    for i in range(len(b)):
        j = np.random.randint(25,NB)
        w[i,j:] = 0.0

    

    t1 = time.perf_counter()
    vols_ref = overlap_volume(torch.tensor(a[0]),torch.tensor(b),torch.tensor(w)).numpy()
    t2 = time.perf_counter()

    print('pytorch:', t2-t1)

    vols_cpu = np.zeros((len(a),len(b)), dtype=NPDTYPE_CPP)

    t1 = time.perf_counter()
    a_cpu = np.ones((a.shape[0], a.shape[1], a.shape[2]+1))
    for i in range(a.shape[0]):
        a_cpu[i,:,:3] = a[i,:,:]

    b_cpu = np.ones((b.shape[0], b.shape[1], b.shape[2]+1))
    for i in range(b.shape[0]):
        b_cpu[i,:,:3] = b[i,:,:]
        b_cpu[i,:,3] = w[i,:]



    overlap_cpp(a_cpu,b_cpu,vols_cpu)
    t2 = time.perf_counter()
    print('cpp', t2-t1)


    #a_gpu = a_cpu.astype(dtype=NPDTYPE_CUDA)
    #b_gpu = b_cpu.astype(dtype=NPDTYPE_CUDA)

    #vols_gpu = np.zeros((len(a),len(b)), dtype=NPDTYPE_CUDA)


    #t1 = time.perf_counter()

    #overlap_cuda(a_gpu,b_gpu,vols_gpu)
    #t2 = time.perf_counter()

    #print('cuda', t2-t1)
    
    #print(vols_cpu)
    #print(vols_gpu)
    #print(vols_ref)
    assert(np.allclose(vols_cpu, vols_ref))
    #assert(np.allclose(vols_gpu, vols_ref))


# test gradient between cpp and pytorch
def test_grad():
    # test the grdaients are correct

    a = np.random.rand(10,3)*10
    b = np.random.rand(15,3)*10
    #b=a
    w = np.ones(15)
    w[13:] = 0.0 # make sure some are masked to imitate the padded array 

    q = np.array([0.9659,0.2588,0,0])


    a = a.astype(dtype=NPDTYPE_CPP)
    b = b.astype(dtype=NPDTYPE_CPP)
    w = w.astype(dtype=NPDTYPE_CPP)
    q = q.astype(dtype=NPDTYPE_CPP)

    RT = quaternion_to_rotation_matrix(torch.tensor(q).unsqueeze(0)).detach().numpy()[0]
    b = np.matmul(b, RT.T)

    #print(a)
    #print(b)

    ta = torch.tensor(a)
    tb = torch.tensor(b).unsqueeze(dim=0)
    tw = torch.tensor(w).unsqueeze(dim=0)

    print(ta.shape, tb.shape, tw.shape)



    qI = torch.tensor([1.0,0.0,0.0,0.0], dtype=TORCHDTYPE, requires_grad = True)

    R = quaternion_to_rotation_matrix(qI)
    t=torch.tensor([5.0,5.0,0.0], requires_grad=True)
    rotb = torch.matmul(tb, R.T) + t 


    #print(a)
    #print(rotb)

    #print(rotb)
    loss = - overlap_volume(ta, rotb, tw)
    loss.backward()
    ref_gradq = qI.grad
    ref_gradt = t.grad

    #print(ref_gradq, ref_gradt)
    ref_grad = torch.cat([ref_gradq, ref_gradt])

    #print(ref_grad.detach().numpy())

    a_cpu = np.ones((a.shape[0], a.shape[1]+1))
    b_cpu = np.ones((b.shape[0], b.shape[1]+1))

    a_cpu[:,:3] = a - t.detach().numpy()
    b_cpu[:,:3] = b 
    b_cpu[:,3] = w

    test_grad_cpu = test_gradient_cpu(a_cpu, b_cpu)
    #print(test_grad_cpu)

    # a = a.astype(dtype=NPDTYPE_CUDA)
    # b = b.astype(dtype=NPDTYPE_CUDA)
    # w = w.astype(dtype=NPDTYPE_CUDA)
    # q = qI.detach().numpy().astype(dtype=NPDTYPE_CUDA)
    # t = t.detach().numpy().astype(dtype=NPDTYPE_CUDA)

    # a_gpu = np.ones((a.shape[0], a.shape[1]+1))
    # b_gpu = np.ones((b.shape[0], b.shape[1]+1))

    # a_gpu[:,:3] = a
    # b_gpu[:,:3] = b 
    # b_gpu[:,3] = w

    # test_grad_gpu = test_gradient_gpu(a_gpu,b_gpu,q,t)
    #print(test_grad_gpu)

    assert(np.allclose(ref_grad.numpy(), test_grad_cpu, rtol=1e-04))
    # assert(np.allclose(ref_grad.numpy(), test_grad_gpu, rtol=1e-04))

    

# def test_omptimize():
#     Nmols= 1000
#     NB = 32
#     a = np.random.rand(1,32,3)*10
#     b = np.random.rand(Nmols,NB,3)*10
#     w = np.ones((Nmols,NB))
#     vols_cpu = np.zeros((len(a), len(b)), dtype=NPDTYPE_CPP)


#     a = a.astype(dtype=NPDTYPE_CPP)
#     b = b.astype(dtype=NPDTYPE_CPP)
#     w = w.astype(dtype=NPDTYPE_CPP)

#     print(b.shape)

#     for i in range(len(b)):
#         j = np.random.randint(25,NB)
#         w[i,j:] = 0.0



#     a_cpu = np.ones((a.shape[0], a.shape[1], a.shape[2]+1))
#     for i in range(a.shape[0]):
#         a_cpu[i,:,:3] = a[i,:,:]

#     b_cpu = np.ones((b.shape[0], b.shape[1], b.shape[2]+1))
#     for i in range(b.shape[0]):
#         b_cpu[i,:,:3] = b[i,:,:]
#         b_cpu[i,:,3] = w[i,:]

#     #print(a_cpu)
#     #print(b_cpu)
        
#     a_gpu = copy.deepcopy(a_cpu).astype(dtype=NPDTYPE_CUDA)
#     b_gpu = copy.deepcopy(b_cpu).astype(dtype=NPDTYPE_CUDA)

#     t1 = time.perf_counter()
#     optimize_overlap_cpu(a_cpu,b_cpu,vols_cpu)
#     t2 = time.perf_counter()

#     print('optim cpu:', t2-t1)
    

#     vols_gpu = np.zeros((len(a),len(b)), dtype=NPDTYPE_CUDA)

#     t1 = time.perf_counter()
#     optimize_overlap_gpu(a_gpu,b_gpu,vols_gpu, True)
#     t2 = time.perf_counter()

#     print('optim gpu:', t2-t1)


#     #print(vols_cpu)
#     #print(vols_gpu)

#     diff = np.fabs((vols_cpu - vols_gpu)/vols_cpu)
#     print(np.max(diff), np.mean(diff))

#     assert(np.max(diff)<0.1)


 


if __name__ == "__main__":
    test_vols()
    test_grad()
    # test_omptimize()
