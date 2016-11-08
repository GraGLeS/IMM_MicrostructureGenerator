The input parameter are provided in the form of three files.

Basic settings
==============

| **NumberOfGrains**
|	Specifies how many grains should be synthezised in the RVE. Number is maintained exactly.
| **NumberOfSubgrains**
|	Specifies how many sub-grains each grain should on average contain.
| **NumberOfPointsPerSubGrain**
|	Specifies the discretization of the domain such that it translates into a cube with a volume equal to the
| 	total volume of spheres with diameter NumberOfPointsPerSubGrain with a total of 
|	NumberOfGrains times (the average) NumberOfSubgrains per grain spheres.

| **MicroGenMode**
|	So far keep 0, this activates the so-called E_VORONOI mode.
| **TextureGEN**
|	Set to 1 to use preference orientation. By default it assigns randomly each grain an orientation from the
|	discrete list of orientations specified in **ReadFromFilename**. Thereafter, it categorizes grain according to
| 	which of the texture components detailed in **AdditionalFilename** are least disoriented to the grain orientation.
|	If none of these texture components is closer than within 10degrees the grain is considered of random orientation.
|	By default this categorization specifies the orientation and stored elastic energy scatter with which the sub-grains
|	inherit orientation and stored elastic energy from their parent grain.

Default grain properties
========================
| **StoredElasticEnergyMax**
|	Specifies the maximum value of the range of stored elastic energy that is assigned uniformly randomly for grains
|	that were characterized to be of random character.
| **StoredElasticEnergyMin**
|	Specifies the minimum value, respectively. Together, the arithmetic average of the two defines the mean of a normal distribution.
| **StoredElasticScatterGrain**
|	Specifies the scattering range, and thus the width (sigma) of a normal distribution according to which the dislocation between grains scatters.
| **StoredElasticScatterSubgrain**
|	Specifies the scattering range, and thus the width (sigma) of this normal distribution according to which dislocation density between the sub-grains inside a specific grain scatters.
| **SubgrainOriScatter**
|	Specifies the only parameter of a Rayleigh distribution (usually referred to as sigma, B in matlab) according to which the orientations
| 	of the sub-grains in default grains are sampled such that the disorientation angle of the grain-local sub-grain population is Rayleigh-distributed.


User-defined macrotexture
=========================
| **ReadFromFilename**
|	Name of a user-defined textfile (see formatting below) which identifies a list of orientations from which to sample the grain's orientations.
| **AdditionalFilename**
|	Name of a user-defined textfile (see formatting below) which identifies texture orientation components into which the grains can be categorized.
|	This enables to assign sub-grains specific properties which differ from grain to grain as explained under **TextureGEN**.
| **CrystalStructure**
|	Currently only fcc implemented. So set to 1.

| **PlotDimension**
|	Set literally according to explanation.

| **ExecuteInParallel**
|	Makes sense to have it activated, i.e. set to 1


Plotting the result as an EBSD map
==================================
| **PlotIPF2DSection**
|	Set to 1 to depict the resulting structure. The following four settings defined **relative** dimensions to draw a subset of the domain.
|	By default an IPFZ map with ND = [0,0,1] is utilized considering antipodal symmetry.
| **PlotWindowXMin**
|	Specifies the left boundary of the domain from which to plot. This limitation is necessary as a structure with millions of sub-grains in 2D would require
|	quickly to draw a picture of several thousand pixel in each dimensions which is seldom required nor practical to be handled
| **PlotWindowXMax**
|	Specifies the right boundary.
| **PlotWindowYMin**
|	Specifies the upper boundary.
| **PlotWindowYMax**
|	Specifies the lower boundary.
|
| By default all plots generate a z = 0 slice.


UserDefinedMacrotexture
=======================

The file contains three single-lined headers which are ignored but must not contain more than one newline each, for instance:

| 1.0\n
| ID x y z vol q0 q1 q2 q3\n
| \n

Thereafter as many lines as desired texture components follow with the following format. The last four specify the quaternion representation.
The data elements are expected to be tab-separated, i.e. by "\t".

integer \t double \t double \t double \t double \t double \t double \t double \t double


Preferential orientation
========================

The file must not contain a header and the format of the textfile is as follows::

  phi1	PHI1	phi2	OriScatter	Stored elastic energy mu	Stored elastic energy sigma	Stored elastic energy internal sigma	SizeScaler	//comments

>
By definition columns become separated by a single tab, i.e. "\t" and a line ends with "\n" The **stored elastic energy mu** and **..sigma**, respectively, is interpreted as the mean and the standard deviation (sigma) of a normal distribution according to which the stored elastic energy of the grains which were categorized into the respective preference orientation is set. Usually the in-grain variation of stored elastic energy is less than between grains and hence the scatter for the sub-grains is defined by **..internal sigma**. In effect, the polycrystal can then have a larger stored elastic energy variation per grain but a lower inside. Similarly, the **oriscatter** defines the preference orientation-specific scatter of a Rayleigh-distribution according to which misorientation is spread. The **sizescaler** defines how much larger or smaller the average sub-grain radius should be in grains of this preference orientation in relation to the global expectation volume. For example if a structure with 100 grains and 1000 sub-grains each at a resolution of 15px is synthesized, and certain grains are given a **sizescaler** argument of 0.25, it means that for these grains on average (1/0.25)^3 more sub-grains will be placed in the same unit grain volume. **The sizescaler is required carefully chosen, i.e. by [0.5, 2] as physically a nanostructured deformed grain next to a several micron sub-grain aggregate is very unlikely and furthermore otherwise the discretization of the finest grain will eventually insufficent (<10px) if the average NumberOfGridPoints is not increased!**