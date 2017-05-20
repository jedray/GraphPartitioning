
## Synopsis 

This project provides an parallel implementation of spectral partitioning for very large sparse graphs with MPI.  
 
## Motivation

Graphs are fun! This project started in [SF2568 Parallel Computations for Large-Scale Problems](https://www.kth.se/social/course/SF2568/) at [KTH Royal Institute of Technology](kth.se).

## Requirements

The implemented algorithm require OpenMPI and ressources for parallel computing. The code was compiled, run and tested on [PDC](https://www.pdc.kth.se/). 
The input graph files must have [METIS](http://glaros.dtc.umn.edu/gkhome/views/metis) format.

## Test 

The compilation command in the makefile uses `mpiicc`. To compile, run the command
          
        $> make compile

The execution command in the makefile uses `mpirun`. To execute, run the command

        $> make 
        

