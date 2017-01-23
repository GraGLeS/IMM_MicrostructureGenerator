# IMM_MicrostructureGenerator

The code is desigend to generate very special input data for GraGLeS. 

It handles 2D and 3D microstructures at once. The main goal of the algorithm is to mimic a microstructrue similar to a recovered one after hotrolling hence, one composed of subgrains. The asumption is, starting with elongated grains after cold rolling, hot rolling will result in a much finer subgrain structure, which is assumed to still correlate with the cold rolled one. Thus we compute another tesselation of these elongated grains in parallel considering different physical properties depending on the original grain orientation, such as 

-stored elastic energy + scatter, 
-subgrain orientation scatter, 
-subgrain stored elastic energy scatter, 
-subgrain size.

The code produces two output files. 
-Container.bin
-Microstructure.uds

Both can be directly reloaded using GraGLeS ( GraGLeS_3D/scripts/VoxelizedParameters.xml ).
