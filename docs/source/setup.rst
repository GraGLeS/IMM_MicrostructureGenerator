How to get and compile?
=======================

| The compilation of the program utilizes open-source libraries as well as standard Linux tools.

1. Check your Linux installation for a working installation of **cmake**, **make**.
2. Check for the successful installation of the following libraries: the **libnuma**, **Eigen**, and **voro++**.
3. If interested in improved memory performance, install for instance Jason Evans' **jemalloc** allocator library.

4. Pull the project from the git repository or the source from which it was provided.
5. Make sure there is a **src** containing the source code and that it includes a **build** folder, and the **CMakeLists.txt** file.
6. Make sure that three input files **parameters.xml**, **PreferenceOrientations.txt**, and **UserDefinedMacrotexture.txt** are within **build**.
7. Open a console, dive into the project and into this **build** folder.
8. Only once when setting up a new computer type **cmake ..**. This inspects your system and generates a customized makefile.
9. Compile the program by typing **make** into the console to utilize this makefile.

| If in between the compilation process unrecoverable errors occur:
| Delete everything in the build folder except for the default input files and start over at step 8.


How to get the additional libraries
===================================



