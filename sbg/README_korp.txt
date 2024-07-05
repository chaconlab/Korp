
######################################################
#  --- KORP v1.22 release - November 26th, 2018 ---  #
#  --- http://chaconlab.org/modeling/korp       ---  #
######################################################


----------------------------------------------
KORP (Knowledge based ORientational Potential)
----------------------------------------------

Pre-compiled 64-bit LINUX binaries (v1.22) are available in the bin/ directory.
We recommend adding binaries path to your PATH environment variable or just
copy them to your favorite "bin" directory. 

Please, choose the appropriate release for your system:

 [Program]        Compiler              Libraries     Linkage     Version
 ------------------------------------------------------------------------
  korpe*          Intel icpc (v16.0.1)      -          static       1.22 
  korpe_gcc       GNU gcc     (v4.6.3)      -          static       1.22   

* Intel compiled binary is the fastest alternative (10-20% faster).


##################
# KORPE TUTORIAL #
##################

The basic instructions to compute the KORP energies with "korpe" tool are 
described here (for both loops and complete proteins). 

The released korp6Dv1.bin energy map contains the information from ALL the
36851 protein chains of our training set (listed in PISCES_id90_r3A.txt file).
PLEASE, BE SURE TO PROVIDE A VALID PATH TO korp6Dv1.bin ENERGY MAP.

The PDB protein chains with resolution better than 3.0 A, R-factor better
than 1.0, and 90% maximum sequence identity were obtained from PISCES server
(http://dunbrack.fccc.edu/Guoli/pisces_download.php) on October, 4th 2018.

More info about KORP is available at http://chaconlab.org/modeling/korp

---------------
 Loop modeling
---------------

The KORP energy computation for loop ensembles is fast and easy with "korpe".


> Single run
------------

For a single case run, just provide a PDB file (filaname.pdb) of the protein 
and a Multi-PDB file with all the loops to be assessed (multi_loops.pdb).

#> korpe filename.pdb --loops multi_loops.pdb --score_file path_to_map/korp6Dv1.bin -o basename

In the output basename.txt file you will find a table with the KORP energies.

For example, in the directory rcd6/ you can run:

#> korpe 3CXM.pdb --loops 3CXM_closed.pdb --score_file ../korp6Dv1.bin -o 3CXM_result

and check the 3FHD_result.txt output file.


> Multiple run 
--------------

For a multiple cases run, use a plain text file with the base names of the 
protein PDBs, i.e. a list of the target PDBs file names without the .pdb
extension (filename_list.txt), and the common suffix of the loops Multi-PDB files 
(_loops.pdb) instead.

#> korpe filename_list.txt --loops _loops.pdb --score_file path_to_map/korp6Dv1.bin -o suffix 

For each XXX entry in the filename_list.txt you will find a file name 
XXX_suffix.txt with the corresponding a table of KORP energies.

For example, in the directory rcd6/ you can run:

#> korpe rcd6_ids.txt --loops _closed.pdb --score_file ../korp6Dv1.bin -o output 


*TIPS
-----

- The presence of one N-terminal and one C-terminal residue (anchors)
in the Multi-PDB file (or the _loops.pdb files) is mandatory.

- Both the protein PDB and the loops Multi-PDB must be numerated 
consistently.

- KORP energy is side-chain independent, it only requires the coordinates
of the N, CA, and C mainchain atoms.

- If benchmarking, add the --rmsd flag and provide a complete PDB (with native
loop) to obtain complete KORP energy statistics in the output_score.txt file.
For example, use these comands for single and multiple runs, respectively, in
the rcd6/ directory:

#> korpe 3CXM.pdb --loops 3CXM_closed.pdb --score_file ../korp6Dv1.bin -o 3CXM_result --rmsd

#> korpe korp6_ids.txt --loops _closed.pdb --score_file ../korp6Dv1.bin -o output --rmsd 


------------------
 Protein modeling
------------------

The KORP energy computation for protein models is fast and easy with "korpe", 


> Single case run
-----------------
For a single case run, just provide a PDB file of the protein and the
corresponding energy map. For example, in directory CASP12DCsel20/ run:

#> korpe T0902D1_s432m2.pdb --score_file ../korp6Dv1.bin

The computed energy will be prompted to screen.


> Multiple run
--------------
For a multiple cases run, just input a plain text file with the base
names of the PDBs, i.e. a list of PDB file names without the .pdb
extension. For example, in directory CASP12DCsel20/ run:

#> korpe casp12dcsel20.txt --score_file ../korp6Dv1.bin -o korp6D
                                           

*TIPS
-----

- Use a correct residue numeration in the PDB, it may be important for
bonding/non-bonding contacts discrimination.

- KORP energy is side-chain independent, it only requires the coordinates 
of the N, CA, and C mainchain atoms.


##################
# KORPE timmings #
##################

Benchmark     Targets  Structures* Time^
                [#]        [#]      [s]
----------------------------------------
rcd6            100      100100      33
rcd8            100      100100      42
rcd10           100      100100      53
rcd12           100      100100      59
casp10all        78       25092     163 
casp11all        58       17507     116
casp12all        32       10056      62
casp10best150    78        9232      67
casp11best150    55        7680      54
casp12best150    28        3745      28
casp10sel20      73        1333      11
casp11sel20      46         914       7
casp12sel20      23         470       5
rosetta41        41        4599      17
itasser56        56       24650      74
3DRobot         200       60200     371
rosetta3DR41     41        4141      17
itasser3DR56     56       22456      83
----------------------------------------
* Native structures included.
^ Computation time obtained using the Intel compiled version of korpe and
the provided all-heavy-atoms benchmarks in an old Intel(R) Core(TM) 
i7-950 (3.07GHz) workstation. Note that korpe would be around two times 
faster if you used just the N, CA, C atom coordinates as input (instead 
of the all-heavy-atoms benchmarks).


----------
REFERENCES
----------

Please, cite our work if our tools or benchmarks become useful for your
research.

> KORP energy and benchmarks
----------------------------
Lopez-Blanco JR and Chacon P (2018). KORP: Knowledge-based 6D
potential for protein structure prediction. (to be published)

> Improved RCD method
---------------------
Lopez-Blanco JR, Canosa-Valls AJ, Li Y, and Chacon P (2016). RCD+: Fast loop
modeling server. NAR (DOI: 10.1093/nar/gkw395).

> Original RCD method
---------------------
Chys P and Chacon P (2013). Random coordinate descent with spinor-matrices and
geometric filters for efficient loop closure. J. Chem. Theory Comput.
9:1821-1829.


-------
CONTACT
-------

Please, feel free to contact with us!
(suggestions or bug reports are welcome)

Jose Ramon Lopez-Blanco (PhD.)
jrlopez@iqfr.csic.es

Pablo Chacon (PhD.)
pablo@chaconlab.org

