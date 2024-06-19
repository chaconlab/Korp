# Korp
Knowledge-base ORientational Potential (KORP) utilizes a 6D joint probability and a minimalist representation to outperform state-of-the-art statistical potentials for protein and loop modeling. The basic instructions for computing KORP energies for loops and proteins are described here. 

Please, cite our work if our tools or benchmarks become useful for your research:

López-Blanco J.R. and Chacón P. (2019) KORP: Knowledge-based 6D potential for fast protein and loop modeling. Bioinformatics 35 (17) 3013–3019. <a href="https://academic.oup.com/bioinformatics/article/35/17/3013/5289323"><img src="https://chaconlab.org/images/publications/pubmed.jpg" alt="" align="top" border="0" /></a> <a href="https://chaconlab.org/PDF/2019_bioinfo.pdf"><img src="https://chaconlab.org/images/publications/acrobaticon4.gif" alt="" border="0" /></a>


## Protein modeling

The KORP energy computation for protein model is fast and easy with "korpe"

### Single case

In this case, just provide a PDB file of the protein and the corresponding energy map. For example, in directory CASP12DCsel20/ run:
```
#> korpe T0902D1_s432m2.pdb --score_file ../korp6Dv1.bin
The computed energy will be prompted to screen
```

### Multiple run

For multiple proteins just input a plain text file with the basenames of the PDBs, i.e. a list of PDB file names without the .pdb extension. As an example, in directory CASP12DCsel20/ you can run:
```
#> korpe casp12dcsel20.txt --score_file ../korp6Dv1.bin -o korp6D
```

Tips & tricks
- Check the path of the energy map file korp6Dv1.bin.
- Use a correct residue numeration in the PDB, it may be important for bonding/non-bonding contact discrimination.
- KORP energy is side-chain independent, it only requires the coordinates of the N, CA, and C atoms.

## Loop modeling

### Single run

For a single case run, just provide a PDB file (filename.pdb) of the protein and a Multi-PDB file with all the loops to be assessed (multi_loops.pdb)
```
#> korpe filename.pdb --loops multi_loops.pdb --score_file path_to_map/korp6Dv1.bin -o basename
```
In the output basename.txt file you will find a table with the KORP energies. For example in the directory rcd6/ you can run:
```
#> korpe 3CXM.pdb --loops 3CXM_closed.pdb --score_file ../korp6Dv1.bin -o 3CXM_result
```
and check the 3CXM_result.txt file.

### Multiple run

For a multiple cases run, use a plain text file with the base names of the protein PDBs, i.e. a list of PDB file names without the .pdb extension (filename_list.txt), and the common suffix of the loops Multi-PDB files (_loops.pdb) instead.
```
#> korpe filename_list.txt --loops _loops.pdb --score_file path_to_map/korp6Dv1.bin -o suffix
```
For each XXX entry in the filename_list.txt you will find a file name XXX_suffix.txt with the corresponding a table with the KORP energies. For example, in the directory rcd6/ you can run:
```
#> korpe rcd6_ids.txt --loops _closed.pdb --score_file ../korp6Dv1.bin -o output
```

Tips & tricks:

- The presence of one N-terminal and one C-terminal residue (anchors) in the Multi-PDB file is mandatory
- Both the protein PDB and the loops Multi-PDB must be numerated consistently
- Check the path of the energy map file korp6Dv1.bin
- KORP energy is side-chain independent, it only requires the coordinates of the N, CA, and C mainchain atoms. You can save reading time by input only N,CA,C pdb files. 
- If benchmarking, add the --rmsd flag and provide a complete PDB (with nativeloop) to obtain complete KORP energy statistics in the output_score.txt file. For example, use these commands for single and multiple runs:
```
#> korpe 1h3n12A559.pdb --loops 1h3n12A559_closed.pdb --score_file ../korp_loop_v101.bin -o output --rmsd
#> korpe rcd6_ids.txt --loops _closed.pdb --score_file ../korp6Dv1.bin -o output --rmsd
```


