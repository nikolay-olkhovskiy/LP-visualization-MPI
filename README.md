# BSF-skeleton
The BSF-skeleton is designed for creating parallel programs in C++ using the MPI library. The scope of the BSF-skeleton is cluster computing systems and iterative numerical algorithms of high computational complexity. The BSF-skeleton completely encapsulates all aspects that are associated with parallelizing a program on a cluster computing system.
The theoretical basis of the BSF-skeleton is the BSF (Bulk Synchronous Farm) model of parallel computations, which allows predicting the scalability boundary  of a parallel algorithm/program at an early stage of its development. The scalability boundary is the maximum number of processor nodes that speedup increases to. A detailed description of the BSF parallel computation model can be found in [1].

The repository includes User manual (https://github.com/leonid-sokolinsky/BSF-skeleton/blob/master/User-manual.pdf) that explains in detail how to quickly create correct and efficient parallel programs for computing clusters using the BSF-skeleton.

References

[1] L.B. Sokolinsky, BSF: a parallel computation model for scalability estimation of iterative numerical algorithms on cluster computing
systems, Journal of Parallel and Distributed Computing, 149 (2021): pp. 193–206. https://doi.org/10.1016/j.jpdc.2020.12.009.
