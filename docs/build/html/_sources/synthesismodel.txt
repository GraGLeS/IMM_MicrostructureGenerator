Model Assumptions
=================

.. * A microstructure volume comprises cubic voxel as the smallest discretization unit.
.. * The deformation microstructure can be idealized in contiguous regions with homogeneous properties (orientation, dislocation density).

* The number of sub-grains per grain is an integer and becomes scaled in such a way that the grain with the volume of **NumberOfSubgrains** times the volume of a sphere with diameter **NumberOfGridpointsPerSubgrain** will have **NumberOfSubgrains** sub-grains. As such, small grains have less sub-grains than large grains. As the structure synthesis results in a deterministic but random realization of a point pattern --- and a Poisson point pattern, in particular --- it is expected that the number of sub-grains is not necessarily **NumberOfGrains** times **NumberOfSubgrainsPerGrain**.

* It is crucial to understand that even though the generator operates determinstically, i.e. it generates the same structure given the same settings, it only does so if the same hardware and the same number of OpenMP threads is utilized. This owes to the fact that internally the generator operates two PRNGs. On the one hand a program local generator for all sequentially executed tasks (with seed -3000) and a thread-local PRNG with a seed 2^31 - *omp_get_thread_num() - 1.






**Basic literature covering the physics of recrystallization microstructure evolution:**
 | Cotterill, P., Mould, P. R.
 | Recrystallization and Grain Growth in Metals
 | Surrey University Press, London, 1976
 
 | Humphreys, F. J., Hatherly, M.
 | Recrystallization and Related Annealing Phenomena
 | Pergamon Press, 2003
 | ISBN: 978-0-08-044164-1
 
 | Gottstein, G.
 | Physical Foundations of Materials Science
 | Springer, Berlin, 2010
 | http://dx.doi.org/doi:10.1007/978-3-662-09291-0
 
 | Gottstein, G., Shvindlerman, L. S.
 | Grain Boundary Migration in Metals: Thermodynamics, Kinetics, Applications
 | CRC Press, Boca Raton, 2010
 | ISBN 9781420054354
 
 | Hallberg, H.
 | Approaches to Modeling of Recrystallization
 | Metals, 2011, 1, 16-48
 | http://dx.doi.org/doi:10.3390/met1010016