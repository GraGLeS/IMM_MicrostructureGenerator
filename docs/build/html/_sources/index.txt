.. IMMs PV-Based Microstructure Generator documentation master file, created by
   sphinx-quickstart on Thu Jul 28 18:32:23 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

| This documents an OpenMP-parallelized, ccNUMA-aware, two-staged Poisson-Voronoi tessellation microstructure generator, 
| which discretizes in 2D and 3D a polycrystal comprising grains and subdividing them into an aggragate of grain-local sub-grains. 
| The latter inherit scattered mean-field properties, like orientation, stored elastic energy from their grains.
| The mean-field properties of the grains, in turn, are sampled from a macrotexture and dislocation density distribution. 
| With these capabilities the program aims to instantiate simple heterogeneous model microstructures that describe an idealized 
| sub-grain structure of a well-recovered high-stacking fault energy cubic alloy.
|
| The source code was developed by Christian Mießen and Markus Kühbach at the Institute for Physical Metallurgy and Metal Physics.

.. .. figure:: ../images/MICROSTRUCTUREGENLogo.jpg
..    :scale: 50%
..    :align: left
..    :target: http://www.imm.rwth-aachen.de/index.php?id=88&L=1
   

1. Getting started
^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2
     
   setup
   synthesismodel
   
   
   
2. Input
^^^^^^^^

| Real-world microstructures are fascinating, complex but as such inspire to search for an as close as possible depiction in a computer model.
| At the same time, though, these polycrystalline aggregates realize a specific compromise to a number of possible microstructures which fulfill
| a variety of geometrical, topological constraints to remain a space-filling and to sample specific property distributions. In our generator we 
| withstain from introducing to many hard to control parameter an instead focus on studying the evolution of well-characterizable structures
| and to study first in how far evolutionary trends of these can be identified and casted into reliable meta models to preferential grain coarsening. 
| Here you find the description of the input parameter.

.. toctree::
   :maxdepth: 2
   
   parmconfig
   
   
   
3. Synthesis
^^^^^^^^^^^^
 
| As with every multi-threaded program, it is necessary to decide on how many resources to operate.
| For the specific case of OpenMP one simply has to set proper values for the environment variable:
|
|	**export OMP_NUM_THREADS=<n>**
|
| With *<n>* identifying the number of threads. For instance, running with two threads requires OMP_NUM_THREADS=2
| Thereafter, the simulation is issued by the following command line call:
|
|	**MicrostructureGenerator**
|
| This loads by default the parameters.xml input file, which specifies the (potentially) user-defined and different
| than default name of two additional files from which to interpret orientations and specific settings for orientation scatter. 

| Preliminary benchmarking of a structure with one million grains in 2D on a AMD 8350FX eight core utilizing four threads 
| took approximately 10 minutes to execute.


4. Output
^^^^^^^^^

.. .. figure:: ../images/GenericResultSingleKinetics.jpg
..    :scale: 40%
..    :align: center

| The generator produces the following output files:
|
| 1. **Container.raw**
| It encodes the voxelized microstructure by unique sub-grain IDs. The IDs are in the range of [1,N]. 
| The container encodes implicitly via unsigned int without specifying endianness explicitly.
| For this reason, on RWTH Aachen University cluster's Intel machines LittleEndianess is expected.
|
| 2. **Microstructure.uds**
| It details the properties of each sub-grain. **Both files in combination only describe the structure.**
|
| 3. **MicrostructureDiagnostics.uds**
| It details in addition to the second file the properties of the grains, such that each sub-grain can be traced back to its parent.
|
| 4. **Microstructure.IPFZ.png**
| The IPFZ colored rendition of a z-slice to the domain.
|
| 5. **Parenthood.bin**
| A binary file which encodes pairs of Subgrain/ParentGrainID to enable the identification 
| of those sub-grains which lay close to grain-grain boundaries during subsequent post-processing steps.


Version history
^^^^^^^^^^^^^^^

| **v1.0** first version, 
| 	OpenMP-parallellized, ccNUMA-aware, two-staged Poisson-Voronoi tessellation microstructure generator,
| 	which discretizes in 2D or 3D a polycrystal comprising grains and subdividing them into an aggragate
| 	of grain-local sub-grains. The latter inherit scattered mean-field properties, like orientation,
| 	stored elastic energy from their grains.


References
^^^^^^^^^^

.. toctree::
  :maxdepth: 2
  
  references
  

Funding
^^^^^^^
 | The authors gratefully acknowledge the support from the DFG in the frame of the Reinhart Koselleck project (GO 335/44-1).
 | Furthermore, we acknowledge the support from the FZ Jülich and RWTH Aachen University for granting us computing time 
 | to develop and verify our methods within the frame of JARAHPC projects.
 
 
 
 
Licence
^^^^^^^

 | The project is licenced under the GNU v2.0.


Questions, contributions
^^^^^^^^^^^^^^^^^^^^^^^^
 | Just let us know or contact *markus.kuehbach@rwth-aachen.de*