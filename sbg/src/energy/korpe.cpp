//============================================================================
// Name        : energy.cpp
// Author      : Mon
// Version     :
// Copyright   : Your copyright notice
// Description : Energy computation tool
//============================================================================

/*###########################################################
# NOTES about "korpe", the program to evaluate KORP energy. #
#############################################################

# KORP + RAMA example (don't forget the weight factors --br and --rf):
#> time korpe korp12_ids.txt --loops _closed.pdb --rmsd -e 56 --score_file ~/korp/Korp6Dv1/korp3A90c16f10_korp12no50_6Dr10n36c8b4yz8_map.bin --bf 1.8 --rama_model 3 --rama_file ~/korp/Korp6Dv1/dunbrack.bin --rf 2 -o ramaKorp_rf2

*/

#include "libenergy/include/korpe.h" // KORP energy
#include "libenergy/include/pd2.h" // PD2 bump energy
#include "libenergy/include/tobi.h" // TOBI energy
#include "libenergy/include/rama.h" // Rama energy
#include "libenergy/include/icosa.h" // ICOSA energy
// 2024 remove rasp2
//#include "librasp2/include/rasp2.h" // RASP side chain repacking
#include "time.h" // Santy's timer (parallelization independent)
#include <string>
#include <cmdl/CmdLine.h>

// #########################################################################
// DEFINITIONS ZONE
// #########################################################################
#define FILENAME 300
#define VERSION "v1.28"
#define NAMESIZE 80  // Filename size for Decoy PDBs

// Required by "quicksort"
#define swap(a,b,t) ((t)=(a),(a)=(b),(b)=(t))


// #########################################################################
// GLOBAL VARIABLES ZONE
// #########################################################################

// Compute RMSD from two sets of 3D points (plain float triplets array)
//  a --> First set of coordinates
//  b --> Second set of coordinates
//  n --> Total number of 3D points
float rmsd_coords(float *a, float *b, int n)
{
	float msd = 0.0; // Mean Squared Deviation
	for(int i=0; i < 3*n; i++)
	{
		msd += powf(a[i]-b[i],2);
	}
	msd /= n; // Compute the mean value
	return sqrtf(msd); // and the Square root
}

// Quick-sort recursive routine to sort everything in between `low' <-> `high'
//  arr   --> Input array of values (float) [v0,v1,v2,...vN-1]
//  index --> Input array of indices (int) [0,1,2,...,N-1]
//  low   --> Low pivot point (for recursion)
//  high  --> High pivot point (for recursion)
void quicksort(float *arr, int *index, int low, int high)
{
	int y; // dummy variable for swap ints
	int i = low; // left
	int j = high; // right
	float z = arr[ index[(low + high) / 2] ]; // pivot

	do {
		while(arr[index[i]] < z) // find member above
			i++;

		while(arr[index[j]] > z) // // find member below
			j--;

		if(i <= j) // Swap two elements
		{
            swap(index[i],index[j],y);
			i++;
			j--;
		}
	} while(i <= j);

	// recurse
	if(low < j)
		quicksort(arr, index, low, j);

	if(i < high)
		quicksort(arr, index, i, high);
}

// Get the maximum value
float getmax(float *a, int n)
{
	float max = a[0];
	for(int i=1; i<n; i++)
		if(a[i] > max)
			max = a[i];
	return max;
}

// Get the minimum value
float getmin(float *a, int n)
{
	float min = a[0];
	for(int i=1; i<n; i++)
		if(a[i] < min)
			min = a[i];
	return min;
}

// Count occurrences within the interval: [start,end[
int numinrange(float *a, int n, float start, float end)
{
	int cont = 0;
	for(int i = 0; i<n; i++)
		if(a[i] >= start && a[i] < end)
			cont++;
	return cont;
}

// Reads all lines from a text file and returns the number of lines read (automatic memory allocation)
//  file       --> File name
//  p_lines    --> Pointer to the lines array (it will be allocated automatically)
//  linelength --> Length of each line (number of characters)
int readTextLines(char *file, char ***p_lines, int linelength = 50) // Reading all rows from textfile
{
	FILE *f;
	char myline[1024]; // For reading text files
	char **lines; // lines list
	int nlines = 0; // Number of lines
	int nscan; // Number of successful format conversions from sscanf
	lines = (char **) malloc( sizeof(char *) );

	if( !(f = fopen(file,"r")) )
	{
		fprintf(stderr,"readTextLines> ERROR opening %s. Forcing exit!\n",file);
		exit(2);
	}

	while( fgets(myline,1024,f) )
	{
		//printf("%s\n",myline);
		if(myline[0] != '#')
		{
			// Memory allocations...
			if( !( lines = (char **) realloc( lines, sizeof(char *) * (nlines+1) ) ) )
			{
				fprintf(stderr,"readTextLines> Memory for text list can't be allocated! Forcing exit!\n\n");
				exit(1);
			}

			lines[nlines] = (char *) malloc( sizeof(char) * linelength );

			nscan = sscanf(myline, "%s", lines[nlines]); // Reading line

			// fprintf(stderr,"line %3d: %s\n",nlines+1,lines[nlines]);

			if(nscan != 1)
			{
				fprintf(stderr,"readTextLines> %d successful format conversions from sscanf()! Forcing exit!\n\n",nscan);
				exit(1);
			}

			nlines++;
		}
	}
	fclose(f);

	*p_lines = lines; // output lines
	return nlines;
}

// Read column "col" (1,2,...,N) from file into "p_data" (automatic memory allocation). Returns the number of elements read.
int readColumn(char *file, int col, float **p_data)
{
	char myline[1024];
	int nscan = 0;
	char fmt[256]="";

	float *data; // Output data
	data = (float *) malloc( sizeof(float) ); // Allocate first data

	// Read input text file with PDB file names to be processed
	FILE *f_file;
	f_file	= fopen(file,"r");		// Read input file
	if(f_file == NULL)				// If problem with opening file...
	{
		fprintf(stderr,"%s> Input file %s could not be opened...\nrcd>\n","readColumn",file);
		exit(1);
	}									// File opening check passed...

	for(int i=0; i<col-1; i++)
		strcat(fmt, " %*s"); // to get Score column from input file
	strcat(fmt, " %f");
	fprintf(stderr,"fmt= %s\n",fmt);

	int cont = 0;
	while( fgets(myline,1024,f_file) )
	{
		//printf("%s\n",myline);
		if(myline[0] != '#')
		{
			// Memory allocations...
			if( !( data = (float *) realloc( data, sizeof(float) * (cont+1) ) ) )
			{
				fprintf(stderr,"%s> Memory for output data can't be allocated! Forcing exit!\n\n","readColumn");
				exit(1);
			}

			nscan = sscanf(myline, fmt, &data[cont]);

			cont++;
		}
	}

	*p_data = data; // output column
	return cont; // output number of elements
}

// Read column "col" (1,2,...,N) from file into "p_data" (automatic memory allocation). Returns the number of elements read.
int readColumn(char *file, int col, char ***p_data)
{
	char myline[1024];
	int nscan = 0;
	char fmt[256]="";

	char **data; // Output data
	data = (char **) malloc( sizeof(char *) ); // Allocate first pointer to data

	// Read input text file with PDB file names to be processed
	FILE *f_file;
	f_file	= fopen(file,"r");		// Read input file
	if(f_file == NULL)				// If problem with opening file...
	{
		fprintf(stderr,"%s> Input file %s could not be opened...\nrcd>\n","readColumn",file);
		exit(1);
	}									// File opening check passed...

	for(int i=0; i<col-1; i++)
		strcat(fmt, " %*s"); // to get Score column from input file
	strcat(fmt, " %s");
	fprintf(stderr,"fmt= %s\n",fmt);

	int cont = 0;
	while( fgets(myline,1024,f_file) )
	{
		//printf("%s\n",myline);
		if(myline[0] != '#')
		{
			// Memory allocations...
			if( !( data = (char **) realloc( data, sizeof(char *) * (cont+1) ) ) ) // reallocate memory for each line of data
			{
				fprintf(stderr,"%s> Memory for output data can't be reallocated! Forcing exit!\n\n","readColumn");
				exit(1);
			}

			if( !(data[cont] = (char *) malloc( sizeof(char) * 100) ) ) // allocate memory for each input data
			{
				fprintf(stderr,"%s> Memory for each input data can't be allocated! Forcing exit!\n\n","readColumn");
				exit(1);
			}

			nscan = sscanf(myline, fmt, &data[cont]); // read data into allocated chunk of memory

			if(nscan != 1)
			{
				fprintf(stderr,"%s> Invalid number of items read, nscan=%d. Forcing exit!\n\n","readColumn",nscan);
				exit(1);
			}

			if( !( data[cont] = (char *) realloc( data, sizeof(char) * strlen(data[cont]) ) ) ) // reallocate memory for each data
			{
				fprintf(stderr,"%s> Memory for each data can't be reallocated! Forcing exit!\n\n","readColumn");
				exit(1);
			}

			cont++;
		}
	}

	*p_data = data; // output column
	return cont; // output number of elements
}

// Compute the average of a float array "data" of "ndata" elements
float average(float *data, int ndata)
{
	double dummy = 0.0;
	for(int i=0; i<ndata; i++)
		dummy += data[i];
	dummy /= ndata;
	return (float)dummy;
}
float average(int *data, int ndata)
{
	double dummy = 0.0;
	for(int i=0; i<ndata; i++)
		dummy += data[i];
	dummy /= ndata;
	return (float)dummy;
}

// Compute the sigma of a float array "data" of "ndata" elements. Optionally, the average can be provided if already known.
float sigma(float *data, int ndata, float avg=0.0)
{
	double dummy = 0.0;
	if(avg == 0.0)
		avg = average(data,ndata);
	for(int i=0; i<ndata; i++)
		dummy += pow(data[i]-avg,2);
	dummy /= ndata;
	return (float)sqrt(dummy);
}
float sigma(int *data, int ndata, float avg=0.0)
{
	double dummy = 0.0;
	if(avg == 0.0)
		avg = average(data,ndata);
	for(int i=0; i<ndata; i++)
		dummy += pow(data[i]-avg,2);
	dummy /= ndata;
	return (float)sqrt(dummy);
}

// Compute the median from "data" array of "ndata" elements
float median(float *data, int ndata)
{
	// Sort array of RMSDs
	int *sorted = (int *) malloc(sizeof(int) * ndata); // Allocate memory for the SORTed indices array for RMSDs (sortr)
	for(int u = 0; u < ndata; u++) // Initialize sorted array of RMSDs
		sorted[u] = u; // Load with the sequential indices
	quicksort(data, sorted, 0, ndata-1); // Quick sort algorithm

	if(ndata == 1)
		return data[0]; // The unique value
	else
		if(ndata % 2 == 0) // even
			return (data[ sorted[ndata/2] ] + data[ sorted[(ndata/2)-1] ]) / 2; // Average of the two median values
		else // odd
			return data[ sorted[ndata/2] ]; // The median value
}

// Compute the Pearson's correlation coefficient applied to a sample ("r" factor)
float corrPearson(float *x, float *y, int n)
{
	float ax = average(x,n);
	float ay = average(y,n);
	float sx = sigma(x,n,ax);
	float sy = sigma(y,n,ay);

	double r = 0.0;
	for(int i=0; i<n; i++)
		r += (x[i]-ax)*(y[i]-ay);
	r /= n;

	return (float) (r/(sx*sy));
}

// Compute the z-score of "n" data elements ("data") wrt some reference value ("ref")
float zscore(float *data, float ref, int n)
{
	float avg = average(data, n);
	float sig = sigma(data, n, avg);
	return (avg - ref) / sig;
}


// Removing Molecules with missing atoms according to selected coarse-grained model
//  Output the number of removed Molecules.
//  p_removed --> (optional) Pointer to the array with the indices of the removed loops
int remove_bad(Macromolecule *mol, Macromolecule **p_molok, int cgmodel, int **p_removed = NULL)
{
	Molecule *molec;
	pdbIter *iter_mol;
	Macromolecule *dummy; // Molecule Dummy
	Macromolecule *molok; // Molecule OK (curated)
	int nbad = 0; // Bad Molecules counter
	int *ibad; // Array indices of bad Molecules
	molok = new Macromolecule();

	//Iterador para recorrer segmentos
	iter_mol = new pdbIter( mol );
	dummy = new Macromolecule(); // construct dummy Macromolecule

	int nmol = iter_mol->num_molecule();
	ibad = (int *) malloc( sizeof(int) * nmol ); // Allocate the maximum possible

	// Iter Molecules
	for ( iter_mol->pos_molecule = 0; !iter_mol->gend_molecule(); iter_mol->next_molecule() )
	{
		molec = (Molecule *) iter_mol->get_molecule(); // Get current molecule
		dummy->add(molec); // Add Molecule ot dummy Macromolecule

		if(dummy->format_residues(false, cgmodel) > 0)
		{
			printf( "%s> Molecule %4d has missing atoms wrt %d CG-model\n", "remove_bad>", iter_mol->pos_molecule, cgmodel);
			ibad[nbad] = iter_mol->pos_molecule; // store indices of bad Molecules
			nbad++; // count bad Molecules
		}
		else
			molok->add(molec); // Add molecule to curated Macromolecule (molok)

		dummy->removeAll(); // Remove the contained elements list leaving the contained elements intact
	}
	delete dummy;

	// Dump bad Molecules indices. Grep by "bad_molecules" tag...
	fprintf(stdout,"Dumping the indices (1,2,...,N) of the %d Molecules with missing atoms. Total number of molecules: %d\nBad_Molecules",nbad,nmol);
	for(int i=0; i<nbad; i++)
		fprintf(stdout," %04d",ibad[i]+1);
	fprintf(stdout,"\n");

	if(p_removed == NULL) // default behavior
		free(ibad); // free indices
	else
		*p_removed = ibad; // output indices

	*p_molok = molok;
	return nbad;
}

// Remove "ni" items with selected "indices" from an "na" elements "array"
//  array   --> Input array
//  na      --> Number of array items
//  indices --> Array of sorted indices of black-listed items (always increasing)
//  ni      --> Number of black-listed items (mainly required for cross-checking purposes)
float *remove_items(float *array, int na, int *indices, int ni=0)
{
	if( na <= ni || na <= 0 || ni < 0 )
	{
		fprintf(stderr, "remove_items> Error, wrong input:  na= %d  ni= %d.  Forcing Exit!\n", na, ni);
		exit(2);
	}

//	float *output = (float *) malloc( sizeof(float) * (na-ni) ); // allocate memory for output items
	float *output = (float *) malloc( sizeof(float) * na ); // allocate the maximum memory that would be required for output items

	int j = 0; // index for output array
	int k = 0; // index of indices array
	int nbad = 0; // Count the bad cases found for cross-checking
	for(int i=0; i<na; i++) // screen array
	{
		if(i == indices[k]) // if current index "i" is found in the black-list
		{
			k++; // do nothing and increment the index "k" of the always-increasing array of indices "indices"
			nbad++; // count black-listed cases found
		}
		else // if not found, dump output to "output" array
		{
			output[j] = array[i];
			j++;
		}
	}
	free(array); // free input array

	if( ni > 0 && nbad != ni )
	{
		fprintf(stderr, "remove_items> Error, only %d items were deleted from the %d black-listed items requested.  Forcing Exit!\n", nbad, ni);
		exit(2);
	}

	output = (float *) realloc( output, sizeof(float) * (na-nbad) ); // trim output array to the necessary lenght
	return output; // output array without deleted items
}


// #########################################################################
// MAIN PROGRAM
// #########################################################################
//
int main(int argc, char *argv[])
{
	bool debug = true; // Set "true" to dump some debug info (SET BY PARSER LATER!)
	bool mon_switch = true; // MON: Debuging stuff... Set true to restore original behavior
	char prog[]="korpe";
	std::string temp;
	char dummy[FILENAME];
	char dummy2[FILENAME];
	char mypdb[FILENAME];
	char mainpdb[FILENAME];
	char mapfile[FILENAME]; // Potential energy map filename
	char ramafile[FILENAME]; // Dunbrack's neighbor-dependent PDFs filename for Ramachandran-based potentials
	char loopsfile[FILENAME];
	char rmsdfile[FILENAME];
	char ligandfile[FILENAME];
	char frodockfile[FILENAME]; // Frodock input data file
	char outbase[FILENAME]; // Output files base name
	char outfile[FILENAME]; // Output results file name
	char gpfile[FILENAME]; // Gnuplot file name
	char bypass_file[30]; // Suffix for energy bypass files
	int emodel = 0; // Potential energy model
	int rama_model = 3; // Ramachandran energy model
	bool loops_switch = false;
	bool outfile_switch = false;
	bool docking_switch = false;
	bool frodock_switch = false; // Frodock input enabled for energy evaluation of protein docking poses
	bool many_ligands = false;
	bool rmsd_switch = false;
	bool rmsdfile_switch = false;
	bool anchors_present = true;
	bool gplot_switch = false; // Set True to enable GNUplot output (output *.gp files)
	bool bypass_switch = false; // Set True to by-pass energy computation and read energies from text file
	bool has_native = false; // Compute native loop energy (enat) if there exist native loop
	bool bf_specified = false;
	bool check_missingatoms = true; // Check missing atoms in "loops-processing mode" (it should be "false" in "bypass mode")
	bool rasp_switch = false; // Set "true" to enable loop side-chains repacking with RASP
	bool rama_switch = false; // Set "true" to activate Rama potential (added as bonding contribution to the selected "emodel" potential)
	bool multiseq_switch = false; // Multiple sequences expected for loops...
	bool loopenv_switch = false; // Consider only loop vs. environment energy in KORP models (-e 5 or 56)
	char chain='*';
	char tobi_type = 2; // 1 or 2 for Tobi-1 (ADP-I) or Tobi-2 (ADP-II) energy models
	double elapsed = 0.0; // Elapsed time variable [s]
	double elapsed2 = 0.0; // Elapsed time variable [s]
	double enat = 0.0; // Native loop energy
	float drmsd = 0.5; // RMSD increment
	float ebins[] = { 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00}; // Energy bin size factors
	float bonding_factor = 0.3; // Bonding energy factor (typically 0.3)
	float rama_factor = 10.0; // Ramachandran energy factor (typically 10)
	int ri,rf;
	int cgmodel = 0; // Atomic model (for formatting)
	int nebins = 7; // Number of energy bins
	int col_bypass = 3; // Column index with the pre-computed energies (only required if "bypass_switch" is True)
	int rmsd_num_atoms_per_res = 3; // Number of atoms per residue considered in RMSD calculations (either 3 or 4)
	timerReal timer,timer2; // Real timer (independent of parallelization)

	// float aax[] = {0.0830,0.0134,0.0586,0.0664,0.0407,0.0734,0.0233,0.0579,0.0571,0.0939,0.0180,0.0428,0.0459,0.0371,0.0507,0.0611,0.0554,0.0715,0.0141,0.0357};

	// COMMAND-LINE PARSER:
	using namespace TCLAP;
	CmdLine cmd(argv[0],"   ", VERSION );
	try {

		printf("EXAMPLES: \n"
			   "  Single-run: %s 1alc.pdb --loops 1alc_loops.pdb --noanchors --rmsd --score_file ~/PROCESS/rcd/loco.score -e 1 --debug --gplot -o myoutput\n"
			   "  Multi-run:  %s seok8_ids.txt --loops _closed.pdb --rmsd --score_file ~/prgs/ORDER_AVE/orave4D.bin -e 5 -o orave\n",prog,prog);

		//
		// Define required arguments no labeled
		//

		UnlabeledValueArg<std::string> Input("PDB","Single run: Input PDB file (Receptor PDB for docking)\n"
				"Multiple run: Plain text file (.txt) with one basename of the PDB and loops per row, e.g. 1tca_1 corresponds to PDB file 1tca_1.pdb and loops file 1tca_1<loops_suffix>","default","pdb");
		cmd.add( Input );

		//
		// Optional arguments
		//

		SwitchArg LoopEnv("","loopenv", "Consider only loop vs. environment energy in KORP models (-e 5 or 56). (default=disabled)", false);
		cmd.add( LoopEnv );

		SwitchArg MultiSeq("","multiseq", "Enable multiple sequences for loops (It does not work yet for SLOW energy models in korp). (default=disabled)", false);
		cmd.add( MultiSeq );

		/* Removed 2024
		SwitchArg RASP("","rasp", "Enable loop side-chains repacking with RASP (default=disabled)", false);
		cmd.add( RASP );
        */
		SwitchArg SkipMissingAtoms("","skip_missingatoms", "Disable missing atoms checks. It may be used with \"bypass mode\". (default=disabled)", false);
		cmd.add( SkipMissingAtoms );

		SwitchArg OnlyScore("","only_score", "Only outputs a single <suffix>_score.txt file for all cases. Use it to minimize output in Multi-run mode. (default=disabled)", false);
		cmd.add( OnlyScore );

		ValueArg<std::string> BypassSuffix("","bypass_file", "Bypass file name, or suffix if multi-run mode enabled. (default=_soap.txt)", false, "_soap.txt", "string");
		cmd.add( BypassSuffix );

		ValueArg<int> Bypass("","bypass", "Column to bypass energy (1,2,...Ncols). It enables energy bypass. (default=disabled)",false,3,"int");
		cmd.add( Bypass );

		SwitchArg GPlot("","gplot", "Enable GNUplot output (output *.gp files). (default=disabled)", false);
		cmd.add( GPlot );

		SwitchArg NoAnchors("","noanchors", "Anchors not present in loops. By default anchors should be found in loops Multi-PDB. (default=disabled)", false);
		cmd.add( NoAnchors );

		SwitchArg RMSD3at("","rmsd3at", "Enable if the atomic model ONLY contains N, CA, and C atoms (No carbonylic O) to set the appropriate number of atoms per residue considered for RMSD calculations (otherwise N, CA, C, and O atoms will be considered) (default=disabled)", false);
		cmd.add( RMSD3at );

		ValueArg<std::string> RmsdPDB("","rmsd_pdb", "Single run: PDB file to compute additional RMSDs\n"
				"Multiple run: PDB files suffix <rmsd_suffix> (e.g. _target.pdb) to compute additional RMSDs",false,"MPDB","string");
		cmd.add( RmsdPDB );

		SwitchArg RMSD("","rmsd", "Dump RMSD (default=disabled)", false);
		cmd.add( RMSD );

		SwitchArg Debug("","debug", "Dump debug info (default=disabled)", false);
		cmd.add( Debug );

        ValueArg<float> RamaFactor("","rf", "Ramachandran energy factor, typically 10 goes fine. (Default=10)",false,10.0,"float");
        cmd.add( RamaFactor );

		ValueArg<std::string> RamaFile("","rama_file", "Neighbor-dependent PDFs must be provided here for any Ramachandran-based potential (our dunbrack.bin file)",false,"rama_file","string");
		cmd.add( RamaFile );

        ValueArg<float> BondingFactor("","bf", "Bonding energy factor, typically 0.3 goes fine. (Default=0.3)",false,0.3,"float");
        cmd.add( BondingFactor );

		ValueArg<std::string> ScoreFile("","score_file", "Potential energy map:"
				"\n\t - ICOSA's plain-text format: <AAi><AAj> <face> <dist> <energy> (for energy model 1)"
//				"\n\t - Dunbrack's neighbor dependent PDFs for our Ramachandran potential (our dunbrack.bin file) (for energy model 6)"
				"\n\t - KORP's energy format (3D, 4D, or 6D, for energy model 5) (default)",false,"score_file","string");
		cmd.add( ScoreFile );

		ValueArg<char> Chain("","chain", "Chain ID (default=first chain of loops Multi-PDB)",false,1,"int");
		cmd.add( Chain );

		ValueArg<std::string> OutFile("o","out", "Output files basename for Single run mode or Output files suffix for Multiple run mode.",false,"txt","string");
		cmd.add( OutFile );

		ValueArg<std::string> Frodock("","frodock", "Frodock's output data file (binnary format) to compute energy of the docking poses."
				"The calculated energies will overwrite the energy registry of said data file into the output file (-o option).",false,"dat","string");
		cmd.add( Frodock );

		ValueArg<std::string> LigandFile("","ligand", "Ligand(s) file in PDB(Multi-PDB) format for docking energy models.",false,"MPDB/txt","string");
		cmd.add( LigandFile );

		ValueArg<std::string> LoopsFile("","loops", "Single run: Loops file in Multi-PDB format.\n"
				"Multiple run: Loops files suffix <loops_suffix>, e.g. _loops.pdb",false,"MPDB","string");
		cmd.add( LoopsFile );

		ValueArg<int> RamaModel("r","rama_model", "Ramachandran energy model: "
				 	 	 	 	 	 	 	 	 "\n\t 0 = Pc/1, "
				 	 	 	 	 	 	 	 	 "\n\t 1 = Pc/Pt (weighted Pt), "
				 	 	 	 	 	 	 	 	 "\n\t 5 = Pc/Pt1, "
                                   				 "\n\t 2 = Plcr/Pt (weighted Pt), "
                                   				 "\n\t 6 = Plcr/Pt1, "
												 "\n\t 3 = Plcr/Pc (default)"
												 "\n\t 4 = Plcr/1",false,3,"int");
		cmd.add( RamaModel );

		ValueArg<int> Emodel("e","energy_model", "Energy model: \n\t 1 = ICOSA protein modeling, "
												 "\n\t  2 = PD2 protein bumps energy, "
												 "\n\t  3  = TOBI-1 protein-protein docking (ADP-I)"
												 "\n\t  4 = TOBI-2 protein-protein docking (ADP-II)"
												 "\n\t  5 = KORP (default)"
												 "\n\t  6 = Rama"
												 "\n\t 56 = KORP + Rama",false,5,"int");
		cmd.add( Emodel );

		// Parse the command line.
		cmd.parse(argc,argv);

		// Load variables from parser
		strcpy(mypdb,((temp=Input.getValue()).c_str()));
		strcpy(mapfile,((temp=ScoreFile.getValue()).c_str()));
		rama_factor = RamaFactor.getValue();
		if(RamaFile.isSet())
		{
			strcpy(ramafile,((temp=RamaFile.getValue()).c_str()));
			rama_switch = true;
		}
		if(LoopsFile.isSet())
		{
			strcpy(loopsfile,((temp=LoopsFile.getValue()).c_str()));
			loops_switch = true;

			if(Chain.isSet())
				chain = Chain.getValue();
		}
		if(RmsdPDB.isSet())
		{
			strcpy(rmsdfile,((temp=RmsdPDB.getValue()).c_str()));
			rmsdfile_switch = true;
		}
		if(OutFile.isSet())
		{
			strcpy(outbase,((temp=OutFile.getValue()).c_str()));
			if(!OnlyScore.isSet())
				outfile_switch = true;
			else
				outfile_switch = false; // Only outputs a single <suffix>_score.txt file for all cases. Only for Multi-run mode. (default=disabled)", false);
		}
		if(LigandFile.isSet())
		{
			strcpy(ligandfile,((temp=LigandFile.getValue()).c_str()));
			docking_switch = true;

			// Select single (.pdb) or multiple (.txt) ligand protocol according to ligand filename extension.
			if( strncmp(".txt",ligandfile+strlen(ligandfile)-4,4) == 0 )
				many_ligands = true;
			else
				many_ligands = false;
		}

		if(Frodock.isSet())
		{
			strcpy(frodockfile,((temp=Frodock.getValue()).c_str()));
			frodock_switch = true;
		}

		rama_model = RamaModel.getValue();
		emodel = Emodel.getValue();
		if(emodel == 6 || emodel == 56)
			rama_switch = true;
		if( rama_switch && !RamaFile.isSet() ) // Rama requested but no Rama file found
		{
			fprintf(stdout,"Error! Please, introduce a valid Ramachandran PDFs file (dunbrack.bin) in --rama_file option.\n");
			exit(1);
		}

		if( (emodel == 3 || emodel == 4) && !docking_switch)
		{
			fprintf(stdout,"Error, in TOBI energy model you must provide structure(s) for the ligand(s). Please, check --ligand option.\n");
			exit(1);
		}

		debug = Debug.isSet(); // Dump debug info
        bf_specified = BondingFactor.isSet();
        bonding_factor = BondingFactor.getValue(); // Bonding energy factor (typically 0.3)
		rmsd_switch = RMSD.isSet();
		anchors_present = !NoAnchors.isSet();
        gplot_switch = GPlot.isSet(); // Enable GNUplot output (output *.gp files)

		strcpy(bypass_file,((temp=BypassSuffix.getValue()).c_str()));
		if(Bypass.isSet())
        {
        	bypass_switch = true; // enable bypass
        	col_bypass = Bypass.getValue(); // set column to bypass energies
			fprintf(stdout,"Energies bypass enabled, they will be taken from %s file(s)\n",bypass_file);
        }

		if(SkipMissingAtoms.isSet())
			check_missingatoms = false; // Check missing atoms in "loops-processing mode" (it should be "false" in "bypass mode")
        /* REMOVED 2024
		if(RASP.isSet())
		{
			rasp_switch = true; // Set "true" to enable loop side-chains repacking with RASP
			emodel = 0; // Energy model must be zero
		}
        */

		multiseq_switch = MultiSeq.isSet(); // Multiple sequences expected for loops...

		loopenv_switch = LoopEnv.isSet(); // Consider only loop vs. environment energy in KORP models (-e 5 or 56)

		// Set the number of atoms per residue considered in RMSD calculations (either 3 or 4)
		if( RMSD3at.isSet() )
			rmsd_num_atoms_per_res = 3; // Number of atoms per residue considered in RMSD calculations (either 3 or 4)
		else
			rmsd_num_atoms_per_res = 4;
	}
	catch ( ArgException& e )
	{
		std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl;
	}
	using namespace std;

	if(debug)
		printf("%s> Welcome to the energy calculation tool for macromolecules!\n",prog);

	// Initialize aminoacids (ADP)
	init_aminoacids(Rosseta,pdb);
//	init_aminoacids();

	double ****loco; // ICOSA potential map
	double icosVertices[NVERT][3]; // icosahedron data
	if(emodel==1)	// ICOSA energy
	{
		//	Allocate dynamic array of dimensions: loco[NAASloco][NAASloco][NDIS][NTRI]
		loco = (double ****) malloc( sizeof(double ***) * NAASloco );
		for(int i=0; i<NAASloco; i++)
		{
			loco[i] = (double ***) malloc( sizeof(double **) * NAASloco );
			for(int j=0; j<NAASloco; j++)
			{
				loco[i][j] = (double **) malloc( sizeof(double *) * NDIS );
				for(int k=0; k<NDIS; k++)
					loco[i][j][k] = (double *) malloc( sizeof(double) * NTRI );
			}
		}

		if(debug)
			printf("%s> Reading ICOSA potential energy map\n",prog);

		read_loco(loco, mapfile); // Load loco energy
		init_icos(icosVertices); // Load an icosaedron of size = 1.0
	}

	// KORP input parameters directly taken from potential map header
	FILE *f_map; // Binary potential map file (it includes Mesh info too)

	korp *map; // KORP's map structure with all stuff required...
	if( emodel == 5 || emodel == 56) // KORP or Rama+KORP
	{
		if(bf_specified) // Overwrite map's bonding_factor with user-defined bonding_factor if parser's bonding_factor has been provided
		{
			// Use parser's bonding_factor
			fprintf(stdout,"%s> Overwrite map bonding_factor with user-defined bonding_factor: %.3f\n",prog,bonding_factor);
			map = readKORP(mapfile,bonding_factor); // Wrapper to read a whole KORP energy map (3/4/6D)
			fprintf(stdout,"%s> Current bonding_factor: %.3f\n",prog,map->bonding_factor);
		}
		else
			map = readKORP(mapfile); // Wrapper to read a whole KORP energy map (3/4/6D)
	}

	// Reading Dunbrack's PDFs
	int size = 72; // Number of bins per dimension of the 2D PDF
	float step = 5; // Angular increment [deg] (may depend on Dunbrack's data set)
	float *****pdfs; // Dunbrack's PDFs
	float *****allmaps; // Complete Ramachandran potential (20x20x20 = 8000 72x72 maps)
	if(emodel == 6 || rama_switch) // Rama potential
	{
		// Load Dunbrack's neighbor-dependent PDFs for Ramachandran-based potentials MAP (Dunbrack_TCB)
		pdfs = read_dunbrack(ramafile,&size,&step); // Reads binary
		// Rama stuff for protein modeling (many single runs)
		fprintf(stdout,"%s> Generating the complete Ramachandran potential (20x20x20 = 8000 72x72 maps)... ",prog);
		fflush(stdout);

		allmaps = gen_allrama(pdfs,rama_model);

		fprintf(stdout,"Done!\n");
		fflush(stdout);
	}

	// Selecting appropriate atomic model (Coarse-Graining model) for selected energy model
	switch(emodel)
	{
	case 1: // for ICOSA
	case 6: // for our Ramachandran potential
		cgmodel = 0; // 0= N,CA,C model
		break;

	case 2: // PD2
		cgmodel = 1; // 1= C5
		break;

	case 3:
		cgmodel = 2; // 2= HA
		break;

	case 4:
		cgmodel = 2; // 2= HA
		break;

	case 5: // KORP
		cgmodel = map->model;
		break;

	case 56: // KORP
		if( map->model != 0)
		{
			fprintf(stderr,"Error, KORP + Rama potential has not been implemented yet for multiple frames or non NCAC atomic models. Forcing exit!\n");
			exit(2);
		}
		cgmodel = 0; // 0= N,CA,C model
		break;
	}

	if(rasp_switch) // Re-packing loop side-chains with RASP
	{
		// cgmodel = 2; // HA, all heavy atoms
		cgmodel = 0; // CA, minimal backbone atoms
	}

	// N,CA,C  Conditions
	Condition *ncac;
	ncac = new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
	ncac->add(" N  ");
	ncac->add(" CA ");
	ncac->add(" C  ");
	Conditions *ncac2= new Conditions();
	ncac2->add(ncac);

	// N,CA,C,O  Conditions (for RMSD calculations)
	Condition *ncaco;
	ncaco = new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
	ncaco->add(" N  ");
	ncaco->add(" CA ");
	ncaco->add(" C  ");
	ncaco->add(" O  ");
	Conditions *ncaco2= new Conditions();
	ncaco2->add(ncaco);

	// N,CA,C,O,CB  Conditions
	Condition *ncacocb;
	ncacocb = new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
	ncacocb->add(" N  ");
	ncacocb->add(" CA ");
	ncacocb->add(" C  ");
	ncacocb->add(" O  ");
	ncacocb->add(" CB ");
	Conditions *ncacocb2= new Conditions();
	ncacocb2->add(ncacocb);

	//
	// Single or Multiple run mode (either processing one PDB or many PDBs at once)
	//
	int npdbs = 1; // Total number of cases to be processed (1 by default)
	char **pdbnames = NULL; // List of PDBs for Multiple Run mode
	char *currentloops; // Current Loops filename
	char **reclist; // Receptor list
	char **liglist; // Ligand list
	char **datlist; // Data list
	if(strncmp(mypdb+strlen(mypdb)-4,".txt",4) == 0)
	{
		fprintf(stdout,"MULTIPLE RUN MODE (%s has .txt extension)\n",mypdb);

		if(frodock_switch) // Frodock stuff
		{
			int ncol1,ncol2,ncol3;
			ncol1 = readColumn(mypdb, 1, &reclist); // Load receptor list
			ncol2 = readColumn(mypdb, 1, &liglist); // Load ligand list
			ncol3 = readColumn(mypdb, 1, &datlist); // Load Frodock's data files list

			// Sanity check for non-self-consistent input
			if( !(ncol1 == ncol2 && ncol2 == ncol3) )
			{
				fprintf(stderr,"Error, check input file. Different number of elements per column found! Forcing exit!\n");
				exit(1);
			}
		}
		else // Protein modeling stuff
		{
			npdbs = readTextLines(mypdb,&pdbnames,NAMESIZE); // Reading all rows from text-file (Multi-run)
			fprintf(stdout,"Read %d cases from %s file\n", npdbs, mypdb);
			//		for(int i=0; i<npdbs; i++)
			//			fprintf(stderr,"Case %3d: %s\n",i+1, pdbnames[i]);
		}
	}
	else if(strncmp(mypdb+strlen(mypdb)-4,".pdb",4) == 0)
	{
		fprintf(stdout,"SINGLE RUN MODE (%s has .pdb extension)\n",mypdb);
		npdbs = 1; // Single-run
		pdbnames = (char **) malloc( sizeof(char *) );
		sprintf(mainpdb,"%s",mypdb);
		pdbnames[0] = mainpdb;
		currentloops = loopsfile;
		sprintf(outfile,"%s.txt",outbase); // Setting single output results file name
		// sprintf(outfile,"%s_%s.txt",pdbnames[0],outbase); // Setting single output results file name
		sprintf(gpfile,"%s.gp",outbase); // Setting single output results file name
	}
	else
	{
		fprintf(stdout,"%s has invalid extension: %s\n Forcing exit!!!\n", mypdb, mypdb+strlen(mypdb)-4);
		exit(2);
	}

	//
	// Loading Loops Multi-PDB
	//
	Macromolecule *loopsr; // Loops as read
	Macromolecule *loops; // Loops according to selected atomic model
	Macromolecule *loopsrmsd; // Loops according to RMSD atomic model
	Macromolecule *environ; // Loop environment
	float *lcoord; // Loop coordinates
	float *coordloopsrmsd; // Loop coordinates for RMSD
	int natoms_loops;
	int nloops;
	int nres_loops; // Total number of residues in loops Multi-PDB
	char *seql; // Loop sequence
	int *iseql; // Loop sequence indices
	char *seqls; // All loops sequence
	int *iseqls; // All loops sequence indices
	int nresl; // Loop number of residues
	int natomsl; // Loop number of atoms

	int nbmax = 10; // Number of N-best cases maximum fraction
	float fb[] = {0.001, 0.002, 0.005, 0.010, 0.020, 0.050, 0.100, 0.200, 0.500, 1.000}; // Fraction Best
	float **all_rmsds; // All RMSDs array (for statistics computation)
	float **all_nbestrmsds; // All N-best-RMSDs array (for statistics computation)
	int **all_sorte; // All Energy-sorted array of indices (for statistics computation)
	float *all_zscores; // All Z-scores array (for statistics computation)
	float *all_corrs; // All cross-correlations array (for statistics computation)
	float *all_bestrmsds; // All best-RMSDs array (for statistics computation)
	float *all_bestranks; // All best-Ranks array (for statistics computation)
	double *all_energies; // All energies for single-run mode (computing energies for a long complete structures list)
	int *all_nloops; // Number of loops for all cases (the number of loops per case may be different)
	all_rmsds = (float **) malloc(sizeof(float *) * npdbs); // memory allocation
	all_sorte = (int **) malloc(sizeof(int *) * npdbs); // memory allocation
	all_zscores = (float *) malloc(sizeof(float) * npdbs); // memory allocation
	all_corrs = (float *) malloc(sizeof(float) * npdbs); // memory allocation
	all_bestrmsds = (float *) malloc(sizeof(float) * npdbs); // memory allocation
	all_bestranks = (float *) malloc(sizeof(float) * npdbs); // memory allocation
	all_nloops = (int *) malloc(sizeof(int) * npdbs); // memory allocation

	// Warning, for the N-best RMSD... the matrix order is reversed: (nbmax x npdbs)
	all_nbestrmsds = (float **) malloc(sizeof(float *) * nbmax); // memory allocation
	for(int x = 0; x < nbmax; x++)
		all_nbestrmsds[x] = (float *) malloc( sizeof(float) * npdbs );

	if(npdbs > 1 && !loops_switch) // Multiple run mode (and no loops)
	{
		all_energies = (double *) malloc( sizeof(double) * npdbs ); // All energies for single-run mode (computing energies for a long complete structures list)

		outfile_switch = false; // Only outputs a single <suffix>_score.txt file for all cases. Only for Multi-run mode. (default=disabled)", false);
	}

	// Main loop
	for(int n = 0; n < npdbs; n++)
	{
		enat = 0.0; // Reset current native energy
		printf( "%s> %-5d Processing %s (of %d) -----------------------------------------\n", prog, n+1, pdbnames[n], npdbs );

		if(npdbs > 1) // Multiple run mode
		{
			sprintf(mainpdb,"%s.pdb",pdbnames[n]);
			sprintf(dummy, "%s%s",pdbnames[n],loopsfile);
			currentloops = dummy;
			sprintf(outfile,"%s_%s.txt",pdbnames[n],outbase); // Setting single output results file name
			sprintf(gpfile,"%s_%s.gp",pdbnames[n],outbase); // Setting single output results file name
		}

		int *resnumsloop; // Loop residue numbers (only mobile)
		char *reschainidsloop; // Loop chain-ids (only mobile)
		int *removed_loops = NULL; // Indices of the removed loops (only for loop mode)
		int nbadloops = 0; // Number of bad loops (they were removed in loop mode)

		if(loops_switch)
		{
			if(debug)
				printf( "%s> Reading Loops Multi-PDB: %s\n", prog, currentloops );
			loopsr = new Macromolecule(currentloops); // Reading loops Multi-PDB into a Macromolecule (each loop will be a molecule)
			loopsr->readPDB(currentloops);
			if(debug)
				printf( "%s> Deleting Hydrogen atoms (if any)...\n", prog );
			loopsr->deleteHYDS();
			if(debug)
				printf( "%s> Deleting Hetero-atoms (if any)...\n", prog );
			loopsr->delete_heteros();
			if(debug)
				printf( "%s> Deleting Water molecules (if any)...\n", prog );
			loopsr->delete_waters();
			//		printf( "%s> Deleting duplicate atoms within residue (if any)...\n", prog );
			//		loopsr->delete_duplicates(); // Remove duplicate atoms

			// Formating residues of all Loops
			if(debug)
				printf( "%s> Formatting Loops Multi-PDB residues order and checking for missing atoms\n", prog );
			if(loopsr->format_residues(false,cgmodel) > 0 && check_missingatoms)
			{
				printf( "%s> Warning, missing atom(s) found in Loops Multi-PDB! (according to %d CG-model)\n", prog, cgmodel );

				// Removing Molecules with missing atoms according to selected coarse-grained model
				Macromolecule *molok;
				nbadloops = remove_bad(loopsr, &molok, cgmodel, &removed_loops);
				loopsr = molok; // MON: Some memory leaks here

				printf( "%s> Removed %d loops with missing atoms according to CG-model %d\n", prog, nbadloops, cgmodel );
			}
			else
				if(debug)
					fprintf(stderr, "%s> No missing atom(s) in Loops Multi-PDB according to CG-model: %d\n", prog, cgmodel );

			if(debug)
			{
				printf( "%s> Written formated Loops Multi-PDB\n", prog );
				loopsr->writePDB( "loopsformatted.pdb" );
			}

			// If not provided, it gets chain ID from first loop
			pdbIter *iter_molec = new pdbIter(loopsr);
			pdbIter *iter_seg;
			Segment *seg;
			Chain *ch;
			if(chain=='*' || npdbs > 1) // If no-chain was specified or in multi-pdb run mode
			{
				iter_molec->pos_chain=0;
				ch=iter_molec->get_chain();
				chain = ch->getName()[0];
			}

			// Get anchor residue indices from first loop in the Multi-PDB file
			iter_molec->pos_segment=0;
			seg = iter_molec->get_segment();
			iter_seg = new pdbIter(seg);
			iter_seg->pos_fragment=0;
			ri = (iter_seg->get_fragment())->getIdNumber(); // get Nt anchor residue number (PDB)
			iter_seg->pos_fragment=iter_seg->num_fragment()-1;
			rf = (iter_seg->get_fragment())->getIdNumber(); // get Ct anchor residue number (PDB), i.e. the last that moves...
			if(!anchors_present) // If anchors are not present in loops, indices must be modified
			{
				ri--;
				rf++;
			}
			delete iter_molec;
			delete iter_seg;
			printf( "%s> Anchor residues obtained from Multi-PDB loop %s --> Nt %d, Ct %d (PDB numeration) Chain_id= \"%c\"\n", prog, currentloops, ri, rf, chain );

			// Mobile loop conditions
			Condition *mobile;
			Conditions *mobile2 = new Conditions();
			if(emodel == 6 || emodel == 56 ) // Rama or KORP+Rama potential (including anchors)
				mobile = new Condition(-1,-1,-1,-1,-1,-1,ri,rf,-1,-1); // get residues from Nt to Ct, i.e.anchors + mobile...
			else // Non-Rama potentials (without anchors)
				mobile = new Condition(-1,-1,-1,-1,-1,-1,ri+1,rf-1,-1,-1); // get residues from Nt+1 to Ct-1, i.e. only those mobile...
			mobile2->add(mobile);

			printf( "%s> Selecting mobile N, CA, C, and O atoms from loops for RMSD calculations\n", prog );


			loops = loopsr->select_cpy(mobile2); // select only mobile residues
			delete loopsr;
			loopsr = loops; // Now "loopsr" has only mobile residues (all atoms read)

			if(rmsd_switch)
			{
				if(debug)
					printf( "%s> Selecting mobile N, CA, C, and O atoms from loops for RMSD calculations\n", prog );
				// MON: the following two lines can be merged into one...
//				loops = loopsr->select_cpy(mobile2); // select only mobile residues
//				loopsrmsd = loops->select_cpy(ncaco2); // select atoms for RMSD from mobile
				loopsrmsd = loopsr->select_cpy(ncaco2); // select atoms for RMSD from mobile
//				delete loops;
				if(debug)
					loopsrmsd->writePDB("loopsrmsd.pdb");

				// Get loop coordinates for RMSD
				if(debug)
					printf("%s> Getting coordinates of all the %d loop atoms for RMSD calculations.\n", prog, loopsrmsd->get_num_atoms());
				loopsrmsd->coordMatrix( &coordloopsrmsd );
				delete loopsrmsd;
			}

			// Selecting the appropriate CG-model for the loops:
			if(rasp_switch) // Re-packing loop side-chains with RASP
			{
				// N,CA,C selection
				if(debug)
					printf( "%s> Selecting N, CA and C loops atoms\n", prog );
				loops = loopsr->select_cpy(ncac2);
			}
			else // Energy computation
			{
				switch(emodel)
				{
				case 1: // for ICOSA
				case 6: // for our Ramachandran potential
					// N,CA,C selection
					if(debug)
						printf( "%s> Selecting N, CA and C loops atoms\n", prog );
					loops = loopsr->select_cpy(ncac2);
					break;

				case 2: // PD2
					// N,CA,C,O,CB selection
					if(debug)
						printf( "%s> Selecting N, CA, C, O and CB loops atoms\n", prog );
					loops = loopsr->select_cpy(ncacocb2);
					break;

				case 5: // KORP
					if(map->model == 0 && mon_switch)
					{    // CA
						if(debug)
							printf( "%s> Selecting N, CA, and C loops atoms\n", prog );
						loops = loopsr->select_cpy(ncac2);
					}
					else if(map->model == 1) // C5
					{
						if(debug)
							printf( "%s> Selecting N, CA, C, O and CB loops atoms\n", prog );
						loops = loopsr->select_cpy(ncacocb2);
					}
					else
					{
						fprintf(stderr,"\nError, invalid CG-model for KORP: %d. Forcing exit!\n\n",map->model);
						exit(2);
					}
					break;

				case 56: // KORP+Rama
					if(map->model == 0 && mon_switch)
					{    // CA
						if(debug)
							printf( "%s> Selecting N, CA, and C loops atoms\n", prog );
						loops = loopsr->select_cpy(ncac2);
					}
					else // MON: this may not be necessary...(above there's something doing the same?)
					{
						fprintf(stderr,"Error, KORP + Rama potential has not been implemented yet for multiple frames or non NCAC atomic models. Forcing exit!\n");
						exit(2);
					}
					break;

				default:
					fprintf(stderr,"\nPlease, introduce a valid CG-model (-e %d is not valid). Forcing exit!\n\n",emodel);
					exit(2);
					break;
				}
				// if(debug)
				//	loops->info(stdout);
			}
			delete loopsr; // read pdb is not needed anymore???? RMSDs????

			//			loopsr = loops->select_cpy(mobile2);
			//			delete loops;
			//			loops = loopsr;
			if(debug)
				loops->writePDB("loopsmobile.pdb");

			iter_molec = new pdbIter(loops);
			nloops = iter_molec->num_molecule(); // each loop is a molecule
			delete iter_molec;

			natoms_loops = loops->get_num_atoms();
			nres_loops = loops->get_num_fragments(); // Total number of residues in loops Multi-PDB

			// First molecule conditions
			Condition *firstmolec;
			firstmolec = new Condition(0,0,-1,-1,-1,-1,-1,-1,-1,-1); // First loop selection
			Conditions *firstmolec2= new Conditions();
			firstmolec2->add(firstmolec);
			Macromolecule *molec0;
			if(debug)
				printf( "%s> Select first loop to get sequence\n", prog );
			molec0 = loops->select_cpy(firstmolec2);
			if(debug)
				molec0->writePDB("single_loop.pdb");

			// Get sequence in 1-letter code
			seql = molec0->get_sequence(); // 1-letter code
			nresl = molec0->get_num_fragments();
			natomsl = molec0->get_num_atoms();
			if(debug)
				printf( "%s> %d residues (%d atoms) loop sequence: (Fasta format)\n> loop\n%s\n", prog, nresl, natomsl, seql );
			if((emodel == 5 || emodel == 56) && map->frame_model == 10 && mon_switch ) // Only required for KORP
			{
				resnumsloop = getResNums(molec0);
				reschainidsloop = getResChainIds(molec0);
//				for(int x=0; x<nresl; x++)
//					fprintf(stderr,"%4d %4d %c\n",x,resnumsloop[x],reschainidsloop[x]);
			}
			delete molec0;

			// Initialize indices sequence (array of indices)
			iseql = (int *) malloc( sizeof(int) * nresl );
			if(multiseq_switch) // Multiple sequences expected for loops...
				iseqls = (int *) malloc( sizeof(int) * nres_loops );

			// Get sequence indices for KORP
			for(int i=0; i<nresl; i++)
				iseql[i] = aa2index(seql[i],aa); // Index for a given AA

			// ... for multiple sequences
			if(multiseq_switch) // Multiple sequences expected for loops...
			{
				seqls = loops->get_sequence(); // 1-letter code

				for(int i=0; i<nres_loops; i++)
					iseqls[i] = aa2index(seqls[i],aa); // Index for a given AA
			}

			if(debug)
			{
				printf( "%s> %d aminoacid Indices: \n", prog, nresl);
				for(int i=0; i<nresl; i++)
					printf( "%d ", iseql[i] );
				printf("\n");
			}

			// Get all loop coordinates
			loops->coordMatrix( &lcoord );
			if(debug)
			{
				printf("%s> Getting coordinates of all the %d loop atoms.\n", prog, natoms_loops);
				loops->writePDB("loops.pdb");
			}

			// Setting cis-Prolines in sequence
			float *dihedrals = (float *) malloc( sizeof(float) * nresl*3 ); // allocate dihedrals memory
			finddihedral( lcoord, nresl*3, dihedrals ); // Get the dihedrals array for supplied coordinates
			if(setCisPro(dihedrals,seql,nresl) > 0)
				fprintf(stdout,"%s> cis-PRO(s) detected: %s\n",prog,seql);
			free(dihedrals); // free dihedrals

			if(debug)
				fprintf(stdout,"%s> End processing %d loops (%c chain) from %s\n",prog,nloops,chain,currentloops);
		}
		// END of loops multi-pdb stuff...

		// Reading main PDB
		Macromolecule *molr = new Macromolecule(mainpdb); // Main PDB
		Macromolecule *mol2,*mol;
		if(debug)
			printf( "%s> Reading PDB: %s\n", prog, mainpdb);
		molr->readPDB(mainpdb); // reading main PDB
		if(debug)
			printf( "%s> Deleting Hydrogen atoms (if any)...\n", prog );
		molr->deleteHYDS();
		if(debug)
			printf( "%s> Deleting Hetero-atoms (if any)...\n", prog );
		molr->delete_heteros();
		if(debug)
			printf( "%s> Deleting Water molecules (if any)...\n", prog );
		molr->delete_waters();

		// Formatting main PDB
		if(debug)
			printf( "%s> Formatting main PDB (receptor) residues order and checking for missing atoms at CG-model %d\n", prog, cgmodel );
		if(molr->format_residues(false,cgmodel) > 0)
		{
			fprintf( stdout, "%s> Error, missing atom(s) found in main PDB (receptor) according to %d CG-model!\n", prog, cgmodel );
			exit(2);
		}
		else
			if(debug)
				fprintf(stderr, "%s> No missing atom(s) in main PDB (receptor)!\n", prog );

		if(debug)
		{
			printf( "%s> Written formated pdb\n", prog );
			molr->writePDB( "formatted.pdb" );
		}

		// Get chain index from main PDB (receptor)
		int ichain=0; // Loop chain index
		pdbIter *iter_molr = new pdbIter(molr);
		Chain *ch;
		for( iter_molr->pos_chain = 0; !iter_molr->gend_chain(); iter_molr->next_chain() ) // screen chains
		{
			ch=iter_molr->get_chain();
			if(ch->getName()[0] == chain)
			{
				ichain = iter_molr->pos_chain; // chain internal index
				break;
			}
		}

		if(debug)
			fprintf(stdout,"Main-PDB chain index for chain %c: ichain= %d\n", chain, ichain);

		Macromolecule *molnatloop; // Native Loop for RMSD calculations
		float *coordnatlooprmsd; // Coordinates of Native Loop for RMSD calculations
		float *coordnatloop; // Coordinates of Native Loop
		mol2 = molr; //

		// Conditions to select native loop for checking
		Condition *natloopsel;
		Conditions *natloopsel2= new Conditions();
		natloopsel = new Condition(-1,-1,-1,-1,ichain,ichain,ri,rf,-1,-1);
		switch(emodel)
		{
		case 1:
		case 6:
		case 56:
			natloopsel->add(" N  ");
			natloopsel->add(" CA ");
			natloopsel->add(" C  ");
			break;
		case 5:
			if(map->model == 0 && mon_switch)
			{   // CA
				natloopsel->add(" N  ");
				natloopsel->add(" CA ");
				natloopsel->add(" C  ");
			}
			else
			{   // C5
				natloopsel->add(" N  ");
				natloopsel->add(" CA ");
				natloopsel->add(" C  ");
				natloopsel->add(" O  ");
				natloopsel->add(" CB ");
			}
			break;
		}
		natloopsel2->add(natloopsel);


		// Checking if native loop is present
		if(loops_switch)
		{
			Macromolecule *natloop; // Native loop
			natloop = molr->select_cpy(natloopsel2,true); // direct selection
			if(natloop->get_num_fragments() != rf-ri+1 && rmsd_switch) // rmsd_switch required because native loop may not be present in real cases.
			{
				fprintf(stderr,"ERROR: Missing residues in Main-PDB's native loop (%d found but %d expected). Forcing exit!\n",natloop->get_num_fragments(),rf-ri+1);
				exit(2);
			}
			int natmissing = 0; // Number of atoms missing in Main-PDB loop
			if( (natmissing = natloop->format_residues(false,cgmodel)) > 0 )
			{
				fprintf(stderr,"ERROR: Missing atoms in Main-PDB's native loop (%d missing). Forcing exit!\n", natmissing);
				exit(2);
			}
			else
				has_native = true; // Native loop exists and is complete for current CG-model

			// Get native loop coordinates (including anchors) for Rama potential
			if(emodel == 6 || emodel == 56)
				natloop->coordMatrix( &coordnatloop );

			delete natloop; // never used again...
		}

		// Remove non-mobile residues from loops
		ri++; // Nt-mobile
		rf--; // Ct-mobile
		if(debug)
			printf( "%s> Select loop residues from Nt_Anchor+1 (ri=%d) to Ct_Anchor-1 (rf=%d) (only mobile)\n", prog, ri, rf );


		// Conditions to select native loop for RMSD calculations
		Condition *loopselrmsd;
		Conditions *loopselrmsd2 = new Conditions();

		// Get Native loop for RMSD computations
		float *coordreflooprmsd; // Coordinates of the reference PDB for additional RMSD computations
		if(rmsd_switch)
		{
			if(debug)
				printf( "%s> Selecting N, CA, C, and O atoms from Native PDB for RMSD calculations (mobile loop)\n", prog );

			loopselrmsd = new Condition(-1,-1,-1,-1,ichain,ichain,ri,rf,-1,-1);
			loopselrmsd->add(" N  ");
			loopselrmsd->add(" CA ");
			loopselrmsd->add(" C  ");
			loopselrmsd->add(" O  ");
			loopselrmsd2->add(loopselrmsd);

			// Getting coordinates single row (pseudo-atom model)
			molnatloop = molr->select_cpy(loopselrmsd2,true); // direct selection
			molnatloop->coordMatrix( &coordnatlooprmsd );
			if(debug)
			{
				printf ("%s> Getting Native Loop coordinates for RMSD computations.\n", prog);
				printf( "%s> Written Native Loop pdb\n", prog );
				molnatloop->writePDB( "natlooprmsd.pdb" );
			}

			// Read external Macromolecule for additional RMSD computations
			if(rmsdfile_switch)
			{
				if(npdbs > 1) // Multiple run mode
					sprintf(dummy, "%s%s",pdbnames[n],rmsdfile); // Join PDB-basename with RMSD-PDBs suffix
				else
					sprintf(dummy, "%s",rmsdfile); // Join PDB-basename with RMSD-PDBs suffix

				fprintf(stderr,"Reading %s\n", dummy);

				Macromolecule *molrmsdr = new Macromolecule(dummy); // PDB for RMSD calculations
				molrmsdr->readPDB(dummy); // reading PDB for RMSD calculations

				// Formatting PDB for RMSDs
				if(debug)
					printf( "%s> Formatting main PDB (receptor) residues order and checking for missing atoms at CG-model %d\n", prog, cgmodel );
				if(molrmsdr->format_residues(false,cgmodel) > 0)
				{
					fprintf( stdout, "%s> Error, missing atom(s) found in PDB for RMSDs according to %d CG-model!\n", prog, cgmodel );
					exit(2);
				}
				else
					if(debug)
						fprintf(stderr, "%s> No missing atom(s) in PDB for RMSDs!\n", prog );

//				fprintf(stderr,"Writting kk.pdb from %s\n", dummy);
//				molrmsdr->writePDB( "kk.pdb" );
//				exit(0);

				Macromolecule *molrmsd  = molrmsdr->select_cpy(loopselrmsd2,true); // direct selection

				if(debug)
					printf( "%s> Reading PDB for RMSD calculations: %s\n", prog, dummy);
				molrmsd->readPDB(dummy); // reading PDB for RMSD calculations
				if(debug)
					printf( "%s> Deleting Hydrogen atoms (if any)...\n", prog );
				molrmsd->deleteHYDS();
				if(debug)
					printf( "%s> Deleting Hetero-atoms (if any)...\n", prog );
				molrmsd->delete_heteros();
				if(debug)
					printf( "%s> Deleting Water molecules (if any)...\n", prog );
				molrmsd->delete_waters();

				molrmsd->coordMatrix( &coordreflooprmsd ); // Get raw array of N, CA, C, and O coordinates

				delete molrmsd; // Not needed anymore
				delete molrmsdr; // Not needed anymore
			}
		}

		// Native Loop selection
		Condition *loopselection;
		Conditions *loopselection2 = new Conditions();
		loopselection = new Condition(-1,-1,-1,-1,ichain,ichain,ri,rf,-1,-1);
		loopselection2->add(loopselection);
		Macromolecule *molnat; // Native macromolecule (complete: environment + loop)

		// Separating Loop from Environment by selection...
		switch(emodel)
		{
//		case 1:
//		case 6:
//			if(has_native)
//			{
//				if(debug)
//					printf( "%s> Selecting N, CA and C atoms for native macromolecule\n", prog );
//				molnat = molr->select_cpy(ncac2);
//				if(debug)
//					printf( "%s> Formatting selected N, CA and C atoms from native macromolecule and checking for missing atoms\n", prog );
//				if(molnat->format_residues(false,cgmodel) > 0)
//					printf( "%s> Warning, missing atom(s) found in native macromolecule!\n", prog );
//			}
//
//			if(loops_switch)
//			{
//				if(debug)
//					printf( "%s> Selecting non-mobile residues in main PDB (ichain= %d), i.e. excluding mobile residues: from %d to %d\n", prog, ichain, ri, rf );
//				mol2 = molr->select_cpy(loopselection2,false); // inverse selection
//				delete molr; // readed pdb is not needed anymore
//			}
//
//			if(debug)
//				printf( "%s> Selecting N, CA and C atoms for environment\n", prog );
//			mol = mol2->select_cpy(ncac2); // WARNING: Native PDB without loop (if loops_switch == true)
//			delete mol2;
//			break;

		case 2: // PD2
			if(debug)
				printf( "%s> Selecting N, CA, C, O and CB atoms\n", prog );
			mol = mol2->select_cpy(ncacocb2); // inverse selection
			delete mol2;
			break;

		case 5: // KORP
		case 1: // ICOSA
		case 6: // Rama
		case 56: // KORP+Rama
			if(emodel == 6 || (map->frame_model == 10 && mon_switch)) // For ICOSA and KORP single frame
			{
				if(has_native)
				{
					if(debug)
						printf( "%s> Selecting N, CA and C atoms for native macromolecule\n", prog );
					molnat = molr->select_cpy(ncac2);
					if(debug)
						printf( "%s> Formatting selected N, CA and C atoms from native macromolecule and checking for missing atoms\n", prog );
//					if(molnat->format_residues(false,0) > 0)
					if(molnat->format_residues(false,cgmodel) > 0)
						printf( "%s> Warning, missing atom(s) found in native macromolecule!\n", prog );
				}

				if(loops_switch)
				{
					printf( "%s> Selecting non-mobile residues in main PDB (ichain= %d), i.e. excluding mobile residues: from %d to %d\n", prog, ichain, ri, rf );
					mol2 = molr->select_cpy(loopselection2,false); // inverse selection
					delete molr; // readed pdb is not needed anymore
				}

				if(debug)
					printf( "%s> Selecting N, CA and C atoms for environment\n", prog );
				mol = mol2->select_cpy(ncac2); // WARNING: Native PDB without loop
			}
//			else
			else if( emodel == 5 )
			{
				if(debug)
					printf( "%s> Selecting N, CA, C, O and CB atoms\n", prog );
				mol = mol2->select_cpy(ncacocb2); // WARNING: Native PDB with loop
			}
			else
			{
				fprintf(stderr,"Error, non-NCAC atomic model not implemented yet for selected potential. Forcing exit!\n");
				exit(2);
			}
			delete mol2; // mol2 is molr
			break;

		case 3:
			tobi_type = 1; // Tobi-1 (ADP-I)
			mol = molr; // Receptor
			break;

		case 4:
			tobi_type = 2; // Tobi-2 (ADP-II)
			mol = molr; // Receptor
			break;

		default:
			mol = molr; // Main PDB
			break;
		}

		int nres = mol->get_num_fragments(); // Number of residues (removed mobile residues may not be considered)
		int natoms = mol->get_num_atoms(); // Number of atoms (removed mobile residues may not be considered)

		int nresnat=0,natomsnat=0; // Native number of residues or atoms

		// Get sequence in 1-letter code
		char *seq = mol->get_sequence(); // PDB sequence in 1-letter code
		if(debug)
			printf( "%s> %d aminoacid sequence: (Fasta format)\n> %s\n%s\n", prog, nres, mainpdb, seq );

		char *seqnat;
		int *iseq, *iseqnat;

//		if(emodel == 1 || (emodel == 5 && map->frame_model == 10 && mon_switch) ) // For ICOSA and KORP single frame
		if(emodel == 1 || ( (emodel == 5 || emodel == 56) && map->frame_model == 10 && mon_switch) ) // For ICOSA and KORP single frame
		{
			// Initialize indices sequence (array of indices)
			iseq = (int *) malloc( sizeof(int) * nres );
			for(int i=0; i<nres; i++)
				iseq[i] = aa2index(seq[i],aa); // Index for a given AA

			if(debug)
			{
				printf( "%s> %d aminoacid Indices: \n", prog, nres);
				for(int i=0; i<nres; i++)
					printf( "%d ", iseq[i] );
				printf("\n");
			}

			if(has_native)
			{
				nresnat = molnat->get_num_fragments(); // Number of residues (removed residues are not considered)
				natomsnat = molnat->get_num_atoms(); // Number of atoms (removed residues are not considered)

				// Get sequence in 1-letter code
				seqnat = molnat->get_sequence(); // 1-letter code

				// Initialize indices sequence (array of indices)
				iseqnat = (int *) malloc( sizeof(int) * nresnat );
				for(int i=0; i<nresnat; i++)
					iseqnat[i] = aa2index(seqnat[i],aa); // Index for a given AA
			}
		}

		// Get Main-PDB coordinates "coord"
		float *coord; // Main-PDB coordinates (loop-less, i.e. environment)
		int *resnums; // Main-PDB residue numbers (loop-less, i.e. environment)
		char *reschainids; // Main-PDB chain-ids (loop-less, i.e. environment)
		float *coordnat; // Main-PDB Native coordinates
		int *resnumsnat; // Main-PDB Native residue numbers
		char *reschainidsnat; // Main-PDB Native chain-ids

//		if(emodel == 1 || emodel == 2 || (emodel == 5 && map->frame_model == 10 && mon_switch) || emodel == 6 || rama_switch) // For ICOSA, PD2, KORP, and Rama
		if(emodel == 1 || emodel == 2 || (emodel == 5 && map->frame_model == 10 && mon_switch) || rama_switch) // For ICOSA, PD2, KORP, and Rama
		{
			// Getting coordinates single row (pseudo-atom model)
			if(debug)
				printf ("%s> Getting main PDB coordinates (at selected coarse-graining model).\n", prog);
			mol->coordMatrix( &coord );

			if(has_native)
				molnat->coordMatrix( &coordnat );
		}

//		if( (emodel == 5 && map->frame_model == 10 && mon_switch) || emodel == 6 || rama_switch ) // Only required for KORP and Rama
		if( (emodel == 5 && map->frame_model == 10 && mon_switch) || rama_switch ) // Only required for KORP and Rama
		{
			resnums = getResNums(mol);
			reschainids = getResChainIds(mol);

			if(has_native)
			{
				resnumsnat = getResNums(molnat);
				reschainidsnat = getResChainIds(molnat);
			}
		}

		// Setting cis-Prolines in sequence only for protein modeling and Rama potential
//		if( !loops_switch && (emodel == 6 || rama_switch) )
		if( !loops_switch && rama_switch )
		{
			float *dihedrals = (float *) malloc( sizeof(float) * nres*3 ); // allocate dihedrals memory
			finddihedral( coord, nres*3, dihedrals ); // Get the dihedrals array for supplied coordinates
			if(setCisPro(dihedrals,seq,nres) > 0)
				fprintf(stdout,"%s> cis-PRO(s) detected: %s\n",prog,seq);
			free(dihedrals); // free dihedrals
		}

		FILE *fout = stdout;
		int anchorCt=0; // C-terminal anchor index (in "mol" Macromolecule)
		int anchorNt=0; // N-terminal anchor index (in "mol" Macromolecule)
		int iCt = -1; // Index (iaa) of Ct anchor residue
		int iNt = -1; // Index (iaa) of Nt anchor residue

		// Coordinates to get N- and C-terminal Ramachandran dihedral angles (required for Rama potentials)
		float coordNt[18]; // Coordinates of the 6 Nt atoms (anchorNt-1 + anchorNt)
		float coordCt[18]; // Coordinates of the 6 Ct atoms (anchorCt + anchorCt+1)

		// Get residue indices of loop ends from main PDB (internal numeration)
		pdbIter *iter_atom;
		Tcoor pos;
		if(loops_switch || rasp_switch)
		{
			pdbIter *iter_chain = new pdbIter(mol);
			pdbIter *iter_frag;
			Fragment *frag;
			char molchain;
			int icont=0; // index counter
			bool exitnow = true;
			for( iter_chain->pos_chain = 0; !iter_chain->gend_chain() && exitnow; iter_chain->next_chain() ) // screen chains
			{
				ch=iter_chain->get_chain();
				molchain = ch->getName()[0];
				iter_frag = new pdbIter(ch);
				if(molchain == chain)
				{
					for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment() && exitnow; iter_frag->next_fragment() ) // screen residues
					{
						frag = iter_frag->get_fragment();
						if(ri-1 == frag->getIdNumber() ) // ri-1 --> Nt anchor (not mobile)
						{
							anchorNt = icont;

							// Get anchorNt atomic coordinates (for Rama potential)
//							if( emodel == 6 || rama_switch )
							if( rama_switch )
							{
								// Get index (iaa) of Nt anchor residue
								for (int j = 0; j < 20;j++)
									if( strncmp(frag->getName(), AA[j].aa_name3, 3) == 0 ) //  seq[i] == AA[j].aa_name3)
									{
										iNt = j; // ID number of Nt anchor
										break;
									}

								// Get anchorNt atomic coordinates (for Rama potential)
								iter_atom = new pdbIter( frag );
								iter_atom->pos_atom = 0;
								( iter_atom->get_atom() )->getPosition(pos); // Get N atom position
								coordNt[ 9] = pos[0];
								coordNt[10] = pos[1];
								coordNt[11] = pos[2];
								iter_atom->pos_atom = 1;
								( iter_atom->get_atom() )->getPosition(pos); // Get CA atom position
								coordNt[12] = pos[0];
								coordNt[13] = pos[1];
								coordNt[14] = pos[2];
								iter_atom->pos_atom = 2;
								( iter_atom->get_atom() )->getPosition(pos); // Get C atom position
								coordNt[15] = pos[0];
								coordNt[16] = pos[1];
								coordNt[17] = pos[2];
								delete iter_atom;
							}

							fprintf(stdout,"%s> Nt anchor residue %s %d (%c chain, iaa=%d) found in main PDB: %d (internal index)\n",prog,frag->getName(),ri-1,chain,iNt,anchorNt);
						}
						if(rf+1 == frag->getIdNumber() ) // rf+1 --> Ct anchor (not mobile)
						{
							anchorCt = icont;

							// Get anchorCt atomic coordinates
//							if( emodel == 6 || rama_switch )
							if( rama_switch )
							{
								// Get index (iaa) of Ct anchor residue
								for (int j = 0; j < 20;j++)
									if( strncmp(frag->getName(), AA[j].aa_name3, 3) == 0 ) //  seq[i] == AA[j].aa_name3)
									{
										iCt = j; // ID number of Ct anchor
										break;
									}

								// Get anchorCt atomic coordinates
								iter_atom = new pdbIter( frag );
								iter_atom->pos_atom = 0;
								( iter_atom->get_atom() )->getPosition(pos); // Get N atom position
								coordCt[0] = pos[0];
								coordCt[1] = pos[1];
								coordCt[2] = pos[2];
								iter_atom->pos_atom = 1;
								( iter_atom->get_atom() )->getPosition(pos); // Get CA atom position
								coordCt[3] = pos[0];
								coordCt[4] = pos[1];
								coordCt[5] = pos[2];
								iter_atom->pos_atom = 2;
								( iter_atom->get_atom() )->getPosition(pos); // Get C atom position
								coordCt[6] = pos[0];
								coordCt[7] = pos[1];
								coordCt[8] = pos[2];
								delete iter_atom;
							}

							fprintf(stdout,"%s> Ct anchor residue %s %d (%c chain, iaa=%d) found in main PDB: %d (internal index)\n",prog,frag->getName(),rf+1,chain,iCt,anchorCt);
							exitnow = false;
						}
						icont++;
					}
				}
				else
					icont += iter_frag->num_fragment();
				delete iter_frag;
			}
			delete iter_chain;
		}

		// Dealing with anchors for loops
//		if( loops_switch && (emodel == 6 || rama_switch) )
		if( loops_switch && rama_switch )
		{
			// Get anchorNt-1 atomic coordinates
			pdbIter *iter_frag = new pdbIter(mol);
			iter_frag->pos_fragment = anchorNt-1;
			iter_atom = new pdbIter( iter_frag->get_fragment() );
			iter_atom->pos_atom = 0;
			( iter_atom->get_atom() )->getPosition(pos); // Get N atom position
			coordNt[0] = pos[0];
			coordNt[1] = pos[1];
			coordNt[2] = pos[2];
			iter_atom->pos_atom = 1;
			( iter_atom->get_atom() )->getPosition(pos); // Get CA atom position
			coordNt[3] = pos[0];
			coordNt[4] = pos[1];
			coordNt[5] = pos[2];
			iter_atom->pos_atom = 2;
			( iter_atom->get_atom() )->getPosition(pos); // Get C atom position
			coordNt[6] = pos[0];
			coordNt[7] = pos[1];
			coordNt[8] = pos[2];
			delete iter_atom;
			// Get anchorCt+1 atomic coordinates
			iter_frag->pos_fragment = anchorCt+1;
			iter_atom = new pdbIter( iter_frag->get_fragment() );
			iter_atom->pos_atom = 0;
			( iter_atom->get_atom() )->getPosition(pos); // Get N atom position
			coordCt[ 9] = pos[0];
			coordCt[10] = pos[1];
			coordCt[11] = pos[2];
			iter_atom->pos_atom = 1;
			( iter_atom->get_atom() )->getPosition(pos); // Get CA atom position
			coordCt[12] = pos[0];
			coordCt[13] = pos[1];
			coordCt[14] = pos[2];
			iter_atom->pos_atom = 2;
			( iter_atom->get_atom() )->getPosition(pos); // Get C atom position
			coordCt[15] = pos[0];
			coordCt[16] = pos[1];
			coordCt[17] = pos[2];
			delete iter_atom;
			delete iter_frag;

			fprintf(stderr,"%s> Anchor Nt-1 and Nt coordinates:\n",prog);
			for(int i=0;i<18;i++)
				fprintf(stderr," %6.2f",coordNt[i]);
			fprintf(stderr,"\n%s> Anchor Ct and Ct+1 coordinates:\n",prog);
			for(int i=0;i<18;i++)
				fprintf(stderr," %6.2f",coordCt[i]);
			fprintf(stderr,"\n");
		}

		// RMSD COMPUTATION
		float *rmsds; // RMSDs array
		rmsds = (float *) malloc(sizeof(float) * nloops); // memory allocation
		if(rmsd_switch) // Computing RMSD wrt Native
		{
			if(debug)
				fprintf(stdout,"Computing RMSDs nresl= %d, rmsd_num_atoms_per_res= %d\n",nresl,rmsd_num_atoms_per_res);

			if( rama_switch ) // Rama potential
			{
				// Rama "coordloopsrmsd" and "coordnatlooprmsd" include both anchors...
				for(int n=0; n<nloops; n++) // Screen loops
				{
					// Anchors must NOT be accounted for in loops
					rmsds[n] = rmsd_coords(coordnatlooprmsd, coordloopsrmsd + (n*nresl+1) * 3 * rmsd_num_atoms_per_res, (nresl-2) * rmsd_num_atoms_per_res); // RMSD(N,CA,C,O)
					// fprintf(stderr,"RMSD %5d   %10f\n",n+1,rmsds[n]);
				}
			}
			else // Non-Rama potentials
			{
				// The "coordloopsrmsd" and "coordnatlooprmsd" in Non-Rama potentials do not include anchors...
				for(int n=0; n<nloops; n++) // Screen loops
				{
					rmsds[n] = rmsd_coords(coordnatlooprmsd, coordloopsrmsd+n*nresl * 3 * rmsd_num_atoms_per_res, nresl * rmsd_num_atoms_per_res); // RMSD(N,CA,C,O)
					// fprintf(stderr,"RMSD %5d   %10f\n",n+1,rmsds[n]);
				}
			}
			free(coordnatlooprmsd); // free nat coords
		}

		float *rmsds2; // RMSDs wrt Reference PDB array
		rmsds2 = (float *) malloc(sizeof(float) * nloops); // memory allocation
		if(rmsdfile_switch) // Computing RMSD wrt Reference PDB
		{
			if(debug)
				fprintf(stdout,"Computing RMSDs wrt Reference PDB, nresl= %d, rmsd_num_atoms_per_res= %d\n",nresl,rmsd_num_atoms_per_res);

			if( rama_switch ) // Rama potential
			{
				// Rama "coordloopsrmsd" and "coordnatlooprmsd" include both anchors...
				for(int n=0; n<nloops; n++) // Screen loops
				{
					// Anchors must NOT be accounted for in loops
					rmsds2[n] = rmsd_coords(coordreflooprmsd, coordloopsrmsd + (n*nresl+1) * 3 * rmsd_num_atoms_per_res, (nresl-2) * rmsd_num_atoms_per_res); // RMSD(N,CA,C,O)
					// fprintf(stderr,"RMSD %5d   %10f\n",n+1,rmsds[n]);
				}
			}
			else // Non-Rama potentials
			{
				// The "coordloopsrmsd" and "coordnatlooprmsd" in Non-Rama potentials do not include anchors...
				for(int n=0; n<nloops; n++) // Screen loops
				{
					rmsds2[n] = rmsd_coords(coordreflooprmsd, coordloopsrmsd+n*nresl*3*rmsd_num_atoms_per_res, nresl*rmsd_num_atoms_per_res); // RMSD(N,CA,C,O)
					// fprintf(stderr,"RMSD %5d   %10f\n",n+1,rmsds[n]);
				}
			}
			free(coordreflooprmsd); // free ref coords
		}

		if(rmsd_switch || rmsdfile_switch) // Current loop coordinates not required anymore
			free(coordloopsrmsd);


		// **************
		// New RASP stuff
		// **************

/*
		if(rasp_switch)
		{
			fprintf(stderr,"RASP START\n");
			// mol->writePDB("mainpdb.pdb");

			// RASP's variables
			CLD *rasp;
			string rasp_mut; // repacking mask (the sequence in 1-char per aminoacid format)
			int rasp_start;
			int rasp_end;
			int rasp_nres;

			rasp_start = anchorNt;
			rasp_end = anchorCt;
			rasp_nres = mol->get_num_fragments();

			fprintf(stderr,"seql: %s\n",seql);

			// Interface between SBG's "Macromolecule" (ours) and RASP's "CLD object" (Zhichao Miao)
			// Converts the "protein" Macromolecule into a new CLD object, inserting dummy residues in the loop region:
			// ("loop_seq") between "start" and "end" residues (PDB indexed)
			rasp = macromol2rasp(mol, rasp_mut, seql, ri-1, rf+1, chain);
			// dump_rasp( rasp );

			fprintf(stderr,"RASP seq: %s\n", rasp_mut.c_str());
			fprintf(stderr,"start= %d  rasp_start= %d  end= %d  rasp_end= %d  rasp_nres= %d\n", ri, rasp_start, rf, rasp_end, rasp_nres);

			// Setting re-packing mask
			for(int i=0; i<=rasp_start; i++)
				rasp_mut[i] = rasp_mut[i] + 32; // this sets to lower-case (upper-case means selected for repacking)
			for(int i=rasp_end; i<rasp_nres; i++)
				rasp_mut[i] = rasp_mut[i] + 32; // this sets to lower-case (upper-case means selected for repacking)
			fprintf(stdout,"Selected loop for RASP: %s\n",rasp_mut.c_str());

			// Setting output Multi-PDB file name
			if(loops_switch) // Multi-PDB file (many structures)
				sprintf(dummy,"%s_%s.pdb",pdbnames[n],outbase);
			else // Single PDB
				sprintf(dummy,"%s.pdb",outbase);

			// Delete previous output file
			FILE *myf = fopen(dummy,"w");
			fclose(myf);

			pdbIter *iter_molec = new pdbIter(loops); // iter loops

			// Get parameters from RASP.ini file
			rasp->GetPar();

			// Set re-packable residues from string
			rasp->GetSub(rasp_mut);

			rasp->GetTop();


				for(int n=0; n<nloops; n++) // Screen loops
				{
					// Interface between SBG's "Macromolecule" (ours) and RASP's "CLD object" (Zhichao Miao)
//					rasp_mut.clear();
//					rasp = macromol2rasp(mol, rasp_mut, seql, ri-1, rf+1);

					// Setting re-packing mask
//					for(int i=0; i<=rasp_start; i++)
//						rasp_mut[i] = rasp_mut[i] + 32; // this sets to lower-case (upper-case means selected for repacking)
//					for(int i=rasp_end; i<rasp_nres; i++)
//						rasp_mut[i] = rasp_mut[i] + 32; // this sets to lower-case (upper-case means selected for repacking)
//					fprintf(stdout,"Selected loop for RASP: %s\n",rasp_mut.c_str());

					// Update the atomic coordinates in RASP's CLD object
					loop2rasp(iter_molec, n, rasp, anchorNt);

					// Rebuild Oxygen atom in residues selected for repacking via "Oss"
					rebuildO(rasp);

					// Populate 4-atoms (N,CA,C,O) backbone "bcbn" (without any checking)
					rasp->Pdb2bcbn();

					// Compute all Phi/Psi dihedral angles from "bcbn"
					rasp->PhiPsi();

					// Load Dunbrack's data directly from file (bbdep11.bin) depending of current Phi and Psi...
					rasp->GetLib(); // Must be here...

//					for(int x=0; x<1000000; x++)
//					{
//					}

					rasp->Build(); // Some initializations? Check whether this can go here...

					rasp->SelfEnergy();
					rasp->PairEnergy();
					rasp->Search();

					append_mdl(dummy, rasp, n+1, rasp_start, rasp_end);

					rasp->bcbn.clear(); // Clear all Backbone
					rasp->sdcn.clear(); // Clear all side-chains
					eraseSC(rasp); // Erase Side-Chains of residues selected for repacking via "Oss" in "Pdb"

					indicator("energy> RASP repacking ", n, nloops);
				}
				indicator("energy> RASP repacking ", nloops, nloops);

			delete rasp;
			delete iter_molec;

			fprintf(stderr,"RASP END\n");

			fprintf(stdout,"Complete loops (with side chains) were dumped into: %s\n",dummy);
			// exit(0);
		}
*/

		//
		// ENERGY COMPUTATION
		//
		Macromolecule *molL;
		double energy = 0.0;
		double envE = 0.0; // Environment vs. Environment energy
		float *energies; // Energies array

		if(bypass_switch) // By-pass energies computation
		{
			if(loops_switch) // Multi-PDB file (many structures)
			{
				//				fprintf(stderr,"n= %d\n", n);
				//				fprintf(stderr,"pdbnames[n]= %s\n", pdbnames[n]);
				//				fprintf(stderr,"bypass_file= %s\n", bypass_file);
				fprintf(stderr,"n= %d  pdbnames[n]= %s   bypass_file= %s\n", n,pdbnames[n], bypass_file);
				if(npdbs > 1)
					sprintf(dummy2,"%s%s",pdbnames[n],bypass_file);
				else
					strcpy(dummy2,bypass_file);
				fprintf(stderr,"By-pass energies from: %s\n",dummy2);

				// Read column "col" from file into "p_data" (automatic memory allocation). It returns the number of elements.
				int nloops2 = readColumn(dummy2, col_bypass, &energies); // Number of energies read from bypass file
				fprintf(stdout,"Read %d energies from %s file\n", nloops2, dummy2);

				// if(removed_loops != NULL)
				if(nbadloops > 0) // If there exist any bad loop...
				{
					fprintf(stderr,"Removed loops indices: ");
					for(int j=0; j<nbadloops; j++)
						fprintf(stderr," %d",removed_loops[j]);
					fprintf(stderr,"\n");

					fprintf(stderr,"nbadloops= %d\n",nbadloops);
					energies = remove_items(energies, nloops2, removed_loops, nbadloops);
					nloops2 -= nbadloops; // update number of loops since the energies of the bad loops were removed from "energies" array

					// Some cross-checking
					if(nloops != nloops2-1)
					{
						// sprintf(dummy, "%s%s",pdbnames[n],loopsfile);
						fprintf(stderr,"Sorry, number of loops mismatch between %s%s (%d) and %s (%d). Forcing exit!\n",pdbnames[n],loopsfile,nloops,dummy2,nloops2);
						exit(2);
					}
					else
						fprintf(stderr,"Number of loops found: %d, and number of energies: %d (the last one should correspond to native energy)\n",nloops,nloops2);
				}


				indicator("energy> By-passing energy calculation",nloops,nloops);

				enat = energies[nloops2-1]; // WATCH OUT! The last one should be the Native energy!!!
				// exit(0);
			}
			else
			{
				fprintf(stderr, "Sorry, senseless options... Forcing exit!\n");
				exit(2);
			}
		}
		else // Standard mode (compute energies)
		{
			if(loops_switch)
				energies = (float *) malloc(sizeof(float) * nloops); // energies memory allocation

			switch(emodel)
			{
			case 1: // ICOSA energy
			{
				if(debug)
					fprintf(stdout,"ICOSA energy model selected!\n");

				if(many_ligands) // If input ligand file is a list of PDB files (.txt)
				{
					if(debug)
						fprintf(stdout,"%s> Multiple ligand protocol selected (single receptor)\n",prog);

					// Read input text file with PDB file names to be processed
					FILE *f_in;
					if( (f_in = fopen(ligandfile,"r")) == NULL)	 // Read input file
					{
						fprintf(stderr,"%s> Input file %s could not be opened... forcing exit!\n", prog, ligandfile);
						exit(1);
					}

					fprintf(stdout, "%s> %7s %-20s %12s %10s\n", prog, "#Ligand", "Name", "ICOSA_energy", "Time [ms]");
					if(outfile_switch)
						fprintf(fout, "%-20s %12s %10s\n","#Ligand_name", "ICOSA_energy", "Time[ms]");

					char myline[1024];
					int nscan = 0;
					int nligands = 0;
					Macromolecule *molL2;

					// strstr(sent, word) != NULL
					while( fgets(myline,1024,f_in) )
					{
						//printf("%s\n",myline);
						if(myline[0] != '#')
						{
							timer.startTimer(); // start timer

							nscan = sscanf(myline,"%s",ligandfile);
							if(nscan != 1)
							{
								fprintf(stderr,"%s> Something went wrong in reading %s ... forcing exit!\n", prog, ligandfile);
								exit(1);
							}
							nligands++;

							if(debug)
								printf( "%s>\n%s> Reading ligand(s) PDB: %s\n", prog, prog, ligandfile );
							molL2 = new Macromolecule(ligandfile);
							molL2->readPDB(ligandfile); // Ligand

							if(debug)
								printf( "%s> Formatting current Ligand %s and checking for missing atoms\n", prog, ligandfile );
							if(int missing = molL2->format_residues(false,cgmodel) > 0)
								printf( "%s> Warning, %d missing atom(s) found in ligand!\n", prog, missing );
							else
								if(debug)
									fprintf(stderr, "%s> No missing atom(s) in ligand!\n", prog );

							if(debug)
								printf( "%s> Selecting N, CA and C atoms\n", prog );
							molL = molL2->select_cpy(ncac2);
							delete molL2; // not needed anymore

							char *seqL;
							int *iseqL;
							int nresL;
							// Get number of ligand residues
							nresL = molL->get_num_fragments();
							// Get sequence in 1-letter code
							seqL = molL->get_sequence(); // 1-letter code
							if(debug)
								printf( "%s> %d aminoacid sequence: (Fasta format)\n> %s\n%s\n", prog, nresL, ligandfile, seqL );

							// Initialize indices sequence (array of indices)
							iseqL = (int *) malloc( sizeof(int) * nresL );
							for(int i=0; i<nresL; i++)
								iseqL[i] = aa2index(seqL[i],aa); // Index for a given AA
							if(debug)
							{
								printf( "%s> %d aminoacid Indices: \n", prog, nresL);
								for(int i=0; i<nresL; i++)
									printf( "%d ", iseqL[i] );
								printf("\n");
							}

							float *coordL;
							// Getting coordinates single row (pseudo-atom model)
							if(debug)
								printf ("%s> Getting ligand PDB coordinates (at selected coarse-graining model).\n", prog);
							molL->coordMatrix( &coordL );

							energy = loco_energy(coord,iseq,nres,coordL,iseqL,nresL,loco,icosVertices); // ICOSA energy calculation

							delete molL; // not needed anymore
							free(coordL);
							free(seqL);
							free(iseqL);

							timer.stopTimer(); // stop timer
							elapsed = timer.getElapsedTime();

							fprintf(stdout, "%s> %7d %-20s %12.6f %10.3f\n", prog, nligands, ligandfile, energy, (float) elapsed * 1000);
							if(outfile_switch)
								fprintf(fout, "%-20s %12.6f %10.3f\n", ligandfile, energy, (float) elapsed * 1000);

						}
					}
					fprintf(stdout,"%s> %d receptor-ligand complexes processed!\n", prog, nligands);
				}
				else // if input ligand file is just a single PDB file
				{

					if(loops_switch) // Multi-PDB file (many structures)
					{
						// Native energy
						enat = loco_energy(coordnat,iseqnat,nresnat,loco,icosVertices);

						// Environment vs. Environment (Constant term)
						envE = loco_energy(coord,iseq,nres,nresl,anchorCt,loco,icosVertices); // Mind the "gap" !!!
						//			fprintf(stderr,"envE= %f\n",envE);

						// Loop vs. Loop  and  Loop vs. Environment
						for(int n=0; n<nloops; n++) // Screen loops
						{
							if(multiseq_switch) // Multiple sequences expected for loops...
								iseql = iseqls + n*nresl; // update indices sequence array

							// The "loco" energy calculation stuff...
							// lindex --> residue index where the loop has been extracted from (in "coord" array), i.e. index of the first right-side residue.
							energies[n] = envE + loco_energy(lcoord+n*nresl*9, iseql, nresl, coord, iseq, nres, anchorCt, loco, icosVertices);

							if(n % 10 == 0)
								indicator("energy> ICOSA energy calculation",n,nloops);
						}
						indicator("energy> ICOSA energy calculation",nloops,nloops);

						delete loops;
						free(lcoord);
					}
					else // Single PDB file with one structure
					{
						energy = loco_energy(coord,iseq,nres,loco,icosVertices);
					}
				}
				break;
			}
			case 2: // PD2 energy
			{
				if(debug)
					fprintf(stdout,"PD2 energy model selected!\n");

				char *type = NULL;
				int *res = NULL;
				int natoms = mol->get_num_atoms(); // number of atoms of the whole PDB

				pd2_type(mol, &type);
				pd2_res(mol, &res);
				if(debug)
					for(int i=0; i<mol->get_num_atoms(); i++)
						fprintf(stderr,"%d ",type[i]);

				if(loops_switch)
				{
					int loopNt,loopCt,loopatoms=0;
					int i=0;
					bool isloop=false;
					while(i<natoms)
					{
						if(res[i]==anchorNt+1 && !isloop)
						{
							loopNt = i;
							isloop = true;
						}
						if(res[i]==anchorCt) // (anchorCt-1)+1 to get the first not-moved atom (of the next residue)
						{
							loopCt = i; // first non-mobile atom of the Ct anchor
							break;
						}
						if(isloop)
							loopatoms++; // Count the total number of atoms in the loop
						i++;
					}
					if(debug)
						fprintf(stderr,"Atomic indices of the loop (%d atoms), Nt= %d and Ct= %d\n",loopatoms,loopNt,loopCt);

					//			// Native (if present) loop energy in the first row of output file
					//			energy = pd2_bump(coord,type,res,natoms);
					//			fprintf(fout,"#Native: %10f\n",energy*100);

					for(int n=0; n<nloops; n++) // Screen loops
					{
						// Updating atomic coordinates of the whole protein
						for(int j = loopNt; j < loopCt; j++)
						{
							coord[j*3]     = lcoord[n*loopatoms*3 + (j-loopNt)*3];
							coord[j*3 + 1] = lcoord[n*loopatoms*3 + (j-loopNt)*3 + 1];
							coord[j*3 + 2] = lcoord[n*loopatoms*3 + (j-loopNt)*3 + 2];
						}

						// The "PD2" energy calculation stuff...
						energies[n] = pd2_bump(coord,type,res,natoms);

						if(n % 10 == 0)
							indicator("energy> PD2 energy calculation",n,nloops);
					}
					indicator("energy> PD2 energy calculation",nloops,nloops);
				}
				else
				{
					energy = pd2_bump(coord,type,res,natoms);
				}

				// Creates an array with the atomic type (N,CA,etc..) for the whole Macromolecule ("PD2" energy calculation stuff)
				// p_type --> Pointer to the "type array" (*p_type == NULL for automatic memory allocation)

				break;
			}
			case 3: // TOBI energy
			case 4:
			{
				if(debug)
					fprintf(stdout,"%s> TOBI-%d energy model selected!\n",prog,tobi_type);

				if(many_ligands) // If input ligand file is a list of PDB files (.txt)
				{
					if(debug)
						fprintf(stdout,"%s> Multiple ligand protocol selected (single receptor)\n",prog);

					// Read input text file with PDB file names to be processed
					FILE *f_in;
					if( (f_in = fopen(ligandfile,"r")) == NULL)	 // Read input file
					{
						fprintf(stderr,"%s> Input file %s could not be opened... forcing exit!\n", prog, ligandfile);
						exit(1);
					}

					int Rnumres, *Rnres, *Rpfirst, *Rcas, *Rtpatom, Lnumres, *Lnres, *Lpfirst, *Lcas, *Ltpatom;
					float *Rcoord, *Lcoord;

					// Getting coordinates (for all atoms and CAs), number of atoms per residue, and atom types for TOBI energy.
					Rnumres = mol->pdbmatrices(&Rcoord, &Rnres, &Rpfirst, &Rcas, &Rtpatom);

					fprintf(stdout, "%s> %7s %-20s TOBI%d_energy %10s\n", prog, "#Ligand", "Name", tobi_type, "Time [ms]");
					if(outfile_switch)
						fprintf(fout, "%-20s TOBI%d_energy %10s\n","#Ligand_name", tobi_type, "Time[ms]");

					char myline[1024];
					int nscan = 0;
					int nligands = 0;

					// strstr(sent, word) != NULL
					timer2.startTimer(); // start timer
					while( fgets(myline,1024,f_in) )
					{
						//printf("%s\n",myline);
						if(myline[0] != '#')
						{
							timer.startTimer(); // start timer

							nscan = sscanf(myline,"%s",ligandfile);
							if(nscan != 1)
							{
								fprintf(stderr,"%s> Something went wrong in reading %s ... forcing exit!\n", prog, ligandfile);
								exit(1);
							}
							nligands++;

							if(debug)
								printf( "%s>\n%s> Reading ligand(s) PDB: %s\n", prog, prog, ligandfile );
							molL = new Macromolecule(ligandfile);
							molL->readPDB(ligandfile); // Ligand

							if(debug)
								printf( "%s> Formatting current Ligand %s and checking for missing atoms\n", prog, ligandfile );
							if(int missing = molL->format_residues(false,cgmodel) > 0)
								printf( "%s> Warning, %d missing atom(s) found in ligand!\n", prog, missing );
							else
								if(debug)
									fprintf(stderr, "%s> No missing atom(s) in ligand!\n", prog );

							// Getting coordinates (for all atoms and CAs), number of atoms per residue, and atom types for TOBI energy.
							Lnumres = molL->pdbmatrices(&Lcoord, &Lnres, &Lpfirst, &Lcas, &Ltpatom);

							// The TOBI energy calculation.
							energy = tobi_energy(Rcoord, Rnumres, Rnres, Rpfirst, Rcas, Rtpatom, Lcoord, Lnumres, Lnres, Lpfirst, Lcas, Ltpatom, tobi_type);

							// Free ligand working arrays
							free(Lnres);
							free(Lpfirst);
							free(Lcas);
							free(Ltpatom);
							free(Lcoord);

							// energy = tobi_energy(mol,molL); // TOBI energy calculation

							timer.stopTimer(); // stop timer
							elapsed = timer.getElapsedTime();

							fprintf(stdout, "%s> %7d %-20s %12.6f %10.3f\n", prog, nligands, ligandfile, energy, (float) elapsed * 1000);
							if(outfile_switch)
								fprintf(fout, "%-20s %12.6f %10.3f\n", ligandfile, energy, (float) elapsed * 1000);

							delete molL; // not needed anymore
						}
					}

					// Free receptor working arrays
					free(Rnres);
					free(Rpfirst);
					free(Rcas);
					free(Rtpatom);
					free(Rcoord);

					timer2.stopTimer(); // stop timer
					elapsed2 = timer2.getElapsedTime();

					fprintf(stdout,"%s> %d receptor-ligand complexes processed in %.2f seconds.\n", prog, nligands, elapsed2);
				}
				else // if input ligand file is just a single PDB file
				{
					if(debug)
					{
						fprintf(stdout,"%s> Single ligand protocol selected (single receptor)\n",prog);
						printf( "%s> Reading ligand PDB: %s\n", prog, ligandfile );
						printf( "%s> Formatting Ligand residues order and checking for missing atoms\n", prog );
					}
					molL = new Macromolecule(ligandfile);
					molL->readPDB(ligandfile); // Ligand

					if(int missing = molL->format_residues(false,cgmodel) > 0)
						printf( "%s> Warning, %d missing atom(s) found in ligand!\n", prog, missing );
					else
						if(debug)
							fprintf(stderr, "%s> No missing atom(s) in ligand!\n", prog );

					energy = tobi_energy(mol, molL, tobi_type); // TOBI energy calculation

					delete molL; // not needed anymore
				}

				break;
			}
			case 5: // KORP energy
			case 56: // KORP+Rama energy
			{
				if(debug)
					fprintf(stdout,"KORP energy model selected!\n");

				int ncont; // Number of contacts
				contact *contacts = NULL;
				pdbIter *iter_molecule;
				Segment *seg;
				Residue *res_loop,*res_pdb;
				Tcoor pos;

				// Compute native loop energy (enat)
				if(has_native)
				{
					// Compute KORP contacts for given PDB structure (automatic memory allocation)
					if(map->frame_model==10 && mon_switch ) // Super fast for single frame model
					{	// Since mol is loop-less, molnat has to be used for native
						frame *framesnat;
						framesnat = frameCoord(coordnat, nresnat, resnumsnat, reschainidsnat);
						ncont = contactCoord(&contacts, map->cutoff, framesnat, nresnat, iseqnat);
						free(framesnat);
					}
					else
					{
						if(map->dimensions == 4)
							ncont = contactPDB_OA(mol, &contacts, map->cutoff, map->iaa, map->mapping, map->nintres);
						else
							ncont = contactPDB(mol, &contacts, map->cutoff, map->iaa, map->mapping, map->nintres);
					}

					// Compute KORP energy
					if(map->dimensions==3) // 3D
						enat = korp3D(contacts,ncont,map);
					else if(map->dimensions==4) // 4D
						enat = korp4D(contacts,ncont,map); // Non-bonding = (nb2 < |i-j|)
					else // 6D
						enat = korp6D(contacts,ncont,map); // Non-bonding = (nb2 < |i-j|)

					free(contacts);
					contacts = NULL; // Enables automatic memory allocation on next iteration

//					fprintf(stderr,"enat= %f  ncont= %d\n",enat,ncont);
					// exit(0);

					// CROSS_TEST...
					//					contact *contacts2 = NULL;
					//					ncont = contactPDB(molnat, &contacts2, cutoff, iaa, mapping, nintres);
					////					molnat->writePDB( "molnat.pdb" );
					//
					//					enat = korp6D(contacts2,ncont,map.maps,meshes,mapping,smapping,fmapping,br,iaa); // Non-bonding = (nb2 < |i-j|)
					//					fprintf(stderr,"enat2= %f\n",enat);
					//					// print_contacts(stderr, contacts2, ncont);
					//					free(contacts2);
					//					contacts2 = NULL; // Enables automatic memory allocation on next iteration

				}

				//				for(int x=0; x<30; x++)
				//				{
				////					if(x%3==0)
				////						fprintf(stderr,"\n");
				//					fprintf(stderr," %10f %10f %10f\n", coord[x], coordnat[x], coord[x] - coordnat[x]);
				//				}
				//				fprintf(stderr,"\n");


				if(loops_switch) // Many loops
				{
					pdbIter *iter_molec = new pdbIter(loops); // iter loops
					pdbIter *iter_pdb = new pdbIter(mol); // iter Main-PDB

					if(map->frame_model==10 && mon_switch) // Super fast for single frame model
					{
						double eenv,eloop,eloopenv;

						// Compute Environment frames
						frame *frames; // Environment frames
						frames = frameCoord(coord, nres, resnums, reschainids);

						// Compute contacts with environment
						ncont = contactCoord(&contacts, map->cutoff, frames, nres, iseq);

						// Compute KORP energy for Environment
						if(map->dimensions==6) // 6D
							eenv = korp6D(contacts,ncont,map); // Non-bonding = (nb2 < |i-j|)
						else if(map->dimensions==3) // 3D
							eenv = korp3D(contacts,ncont,map);
						else // 4D
						{
							fprintf(stderr,"Sorry, 4D KORP not yet implemented! Forcing exit!\n");
							exit(2);
						}
						free(contacts);
						contacts = NULL; // Enables automatic memory allocation on next iteration

//						// OpenMP es una puta maravilla!!! "one-liner 5x speed increase on 8 cores (2011)"
//						#pragma omp parallel for
//						#pragma omp parallel for
						for(int n=0; n<nloops; n++) // Screen loops
						{
							if(multiseq_switch) // Multiple sequences expected for loops...
								iseql = iseqls + n*nresl; // update indices sequence array

							// Computing loop frames
							frame *framesloop; // Loop frames
//							fprintf(stderr,"nresl= %d  \n",nresl);
							if(emodel == 5) // KORP-only energy model has two residues less than Rama-based models
								framesloop = frameCoord(lcoord+n*nresl*9, nresl, resnumsloop, reschainidsloop);
							else
								framesloop = frameCoord(lcoord+(n*nresl+1)*9, nresl-2, resnumsloop+1, reschainidsloop+1);

							// Compute contacts of Loop with itself
							if(emodel == 5) // KORP-only energy model has two residues less than Rama-based models
								ncont = contactCoord(&contacts, map->cutoff, framesloop, nresl, iseql);
							else
								ncont = contactCoord(&contacts, map->cutoff, framesloop, nresl-2, iseql+1);

							// Compute KORP energy for Loop vs. itself
							if(map->dimensions==6) // 6D
								eloop = korp6D(contacts,ncont,map); // Non-bonding = (nb2 < |i-j|)
							else if(map->dimensions==3) // 3D
								eloop = korp3D(contacts,ncont,map);
							else // 4D
							{
								fprintf(stderr,"Sorry, 4D KORP not yet implemented! Forcing exit!\n");
								exit(2);
							}

							free(contacts);
							contacts = NULL; // Enables automatic memory allocation on next iteration

							// Compute contacts for Loop vs. Environment
							if(emodel == 5) // KORP-only energy model has two residues less than Rama-based models
								ncont = contactCoord(&contacts, map->cutoff, frames, nres, iseq, anchorNt, framesloop, nresl, iseql);
							else
								ncont = contactCoord(&contacts, map->cutoff, frames, nres, iseq, anchorNt, framesloop, nresl-2, iseql+1);
							free(framesloop); // Free loop frames

							// Compute KORP energy for Loop vs. Environment
							if(map->dimensions==6) // 6D
								eloopenv = korp6D(contacts,ncont,map); // Non-bonding = (nb2 < |i-j|)
							else if(map->dimensions==3) // 3D
								eloopenv = korp3D(contacts,ncont,map);
							else // 4D
							{
								fprintf(stderr,"Sorry, 4D KORP not yet implemented! Forcing exit!\n");
								exit(2);
							}

							if(loopenv_switch)
								energies[n] = eloopenv; // MON: TESTING!!!
							else
								energies[n] = eloopenv + eloop + eenv;


							//fprintf(stderr,"KORP_ENERGY: %f  %f %f %f\n",energies[n], eloopenv, eloop , eenv);

							free(contacts);
							contacts = NULL; // Enables automatic memory allocation on next iteration

							if(n % 10 == 0)
								indicator("energy> KORP energy calculation",n,nloops);
						}

						free(frames); // Free environment frames
					}
					else // slower version... but valid for all KORP frame-models
					{
						// MON: CHECK THIS FOR KORP+RAMA POTENTIAL...
						for(int n=0; n<nloops; n++) // Screen loops
						{
							// Updating atomic coordinates of the loop within the complete protein (Setting loop coordinates into Main-PDB)
							iter_molec->pos_molecule = n; // current segment is the current loop
							iter_molecule = new pdbIter( (Molecule *) iter_molec->get_molecule() );
							iter_pdb->pos_fragment = anchorNt+1; // Nt-anchor residue + 1 (the first mobile residue of the loop)

							// Screen loop residues
							for( iter_molecule->pos_fragment=0; !iter_molecule->gend_fragment(); iter_molecule->next_fragment() )
							{
								res_loop = (Residue *) iter_molecule->get_fragment();
								res_pdb = (Residue *) iter_pdb->get_fragment();
								fprintf(stderr,"Overwriting coordinates of Main-PDB residue %d with Loop residue %d\n",res_loop->getIdNumber(),res_pdb->getIdNumber());

								// Some checking
								if(res_loop->getIdNumber() != res_pdb->getIdNumber())
								{
									fprintf(stderr,"ERROR: Main-PDB (%d) and Loop (%d) residue indices mismatch! Forcing exit!\n",res_loop->getIdNumber(),res_pdb->getIdNumber());
									exit(1);
								}

								pdbIter *iter_atom_loop = new pdbIter( res_loop );
								pdbIter *iter_atom_pdb = new pdbIter( res_pdb );

								// Screen atoms of current loop residue
								for( iter_atom_loop->pos_atom = 0; !iter_atom_loop->gend_atom(); iter_atom_loop->next_atom() )
								{
									iter_atom_pdb->pos_atom = iter_atom_loop->pos_atom; // select the same atom in the pdb
									(iter_atom_loop->get_atom())->getPosition(pos); // Get loop atom position
									(iter_atom_pdb->get_atom())->setPosition(pos); // Set pdb atom position
								}

								iter_pdb->next_fragment(); // next main-pdb residue
								delete iter_atom_loop;
							}
							delete iter_molecule;
							// mol->writePDB( "loopupdated.pdb" );
							// exit(0);

							// Compute KORP contacts for given PDB structure (automatic memory allocation)
							if(map->dimensions == 4)
								ncont = contactPDB_OA(mol, &contacts, map->cutoff, map->iaa, map->mapping, map->nintres);
							else
								ncont = contactPDB(mol, &contacts, map->cutoff, map->iaa, map->mapping, map->nintres);

							// Compute KORP energy
							if(map->dimensions==3) // 3D
								energies[n] = korp3D(contacts,ncont,map);
							else if(map->dimensions==4) // 4D
								energies[n] = korp4D(contacts,ncont,map); // Non-bonding = (nb2 < |i-j|)
							else // 6D
								energies[n] = korp6D(contacts,ncont,map); // Non-bonding = (nb2 < |i-j|)

							free(contacts);
							contacts = NULL; // Enables automatic memory allocation on next iteration

							if(n % 10 == 0)
								indicator("energy> KORP energy calculation",n,nloops);
						}
					}
					indicator("energy> KORP energy calculation",nloops,nloops);

					delete iter_molec;
					delete iter_pdb;
					delete loops;
					if(emodel != 56) // In "56" mode, lcoord is required later by Rama stuff...
						free(lcoord);
				}
				else
				{
//					mol->writePDB( "kk2.pdb" );
//					exit(0);

					// Compute KORP contacts for given PDB structure (automatic memory allocation)
					if(map->dimensions == 4)
						ncont = contactPDB_OA(mol, &contacts, map->cutoff, map->iaa, map->mapping, map->nintres);
					else
						ncont = contactPDB(mol, &contacts, map->cutoff, map->iaa, map->mapping, map->nintres);

					// fprintf(stderr,"Number of contacts: %d\n",ncont);
					// print_contacts(stderr, contacts, ncont, true);

					// Compute KORP energy
					if(map->dimensions==3) // 3D
						energy = korp3D(contacts,ncont,map);
					else if(map->dimensions==4) // 4D
						energy = korp4D(contacts,ncont,map); // Non-bonding = (nb2 < |i-j|)
					else // 6D
						energy = korp6D(contacts,ncont,map); // Non-bonding = (nb2 < |i-j|)

					free(contacts);
					contacts = NULL; // Enables automatic memory allocation on next iteration
				}

				break;
			}
			case 6: // Rama
			{
				break;
			}
			default:
				if(!rasp_switch)
				{
					fprintf(stdout,"Please, select a valid Energy model! emodel= %d\n",emodel);
					exit(0);
				}

				if(loops_switch)
					for(int n=0; n<nloops; n++) // Screen loops
						energies[n] = 0.0; // Zero energy by default
				else
					energy = 0.0; // Zero energy by default
				break;
			}
		}
		delete mol;

//		if(emodel == 6) // Rama potential (must be here to be combined with other potentials)
		if(rama_switch) // Rama potential (must be here to be combined with other potentials)
		{
			if(debug)
				fprintf(stdout,"Ramachandran energy model selected!\n");

			if(loops_switch) // Many loops
			{
				if(emodel == 6) // Rama potential only
					for(int n=0; n<nloops; n++) // Screen loops
						energies[n] = 0.0; // Zero energy by default

				int *iaa = NULL; // AA index for each loop aminoacid
				iaa = seq2iaa(seql,nresl); // Here, "nresl" is the number of loop residues (including anchors)
				// fprintf(stdout,"seql= %s\n",seql);

				float ***maps=NULL; // PDF-maps for each loop aminoacid
				if( !(maps = (float ***) malloc(sizeof(float **) * nresl) ) )
				{
					printf("Sorry, unable to allocate maps memory!!!\n");
					exit(1);
				}

				//  Assigning neighbor dependent distributions for current loop
				int ileft, iright; // left/right 1-letter code aminoacids
				for(int i=1; i < nresl-1; i++) // Screen LOOP AAs (mobile residues only)
				{
					// fprintf(stderr,"AA=%d  iaa[i]=%d iaa[i-1]=%d  iaa[i+1]=%d  size=%d\n",i,iaa[i],iaa[i-1],iaa[i+1],size);
					if(iaa[i-1]==20) // If cis-Pro to the left... (all Prolines are merged when Left or Right side, see Dunbracks paper...)
						ileft = 12; // left-side cis-Pro is set back to trans-Pro (code 12)
					else
						ileft = iaa[i-1];
					if	(iaa[i+1]==20) // If cis-Pro to the right... (all Prolines are merged when Left or Right side, see Dunbracks paper...)
						iright = 12; // right-side cis-Pro is set back to trans-Pro (code 12)
					else
						iright = iaa[i+1];

					maps[i] = allmaps[iaa[i]][ileft][iright];

					// char kk[100];
					// sprintf(kk,"mierda_%02d.txt",i);
					// show_map(maps[i],size,kk); // Gnuplot command:  plot "mierda_01.txt" using 1:2:3 with image
				}

				// Anchor aminoacids should be treated separately
				// Nt:
				if	(iaa[1]==20) // If cis-Pro to the right... (all Prolines are merged when Left or Right side, see Dunbracks paper...)
					iright = 12; // right-side cis-Pro is set back to trans-Pro (code 12)
				else
					iright = iaa[1];
				maps[0] = allmaps[iaa[0]][iNt][iright];
				// Ct:
				if(iaa[nresl-2]==20) // If cis-Pro to the left... (all Prolines are merged when Left or Right side, see Dunbracks paper...)
					ileft = 12; // left-side cis-Pro is set back to trans-Pro (code 12)
				else
					ileft = iaa[nresl-2];
				maps[nresl-1] = allmaps[iaa[nresl-1]][ileft][iCt];

				free(iaa);
				// exit(0);

				// Allocate memory for current target dihedrals
				float *dihedrals = (float *) malloc( sizeof(float) * nresl*3 );
				float *dihedralsXt = (float *) malloc( sizeof(float) * 9 );

				// Compute native loop energy (enat)
				if(has_native)
				{
					double enatRama;

					// Get the dihedrals array (dihedral) for supplied coordinates (co)
					finddihedral( coordnatloop, nresl*3, dihedrals );

					fprintf(stderr,"Native dihedrals: ");
					for(int x=0; x<nresl*3; x++)
						fprintf(stderr," %6.2f",dihedrals[x]*180/M_PI);
					fprintf(stderr,"\n");

					// Compute KORP energy
//					enat = rama_energy2(dihedrals, nresl, maps, size);
					enatRama = rama_energy2(dihedrals, nresl, maps, size);

					// Nt anchor contribution
					finddihedral( coordNt, 6, dihedralsXt );
					dihedralsXt[6] = dihedrals[3]; // Load the Nt anchor's Psi (mobile)

					fprintf(stdout,"Nt dihedrals: ");
					for(int i=0; i<9; i++)
						fprintf(stdout," %8.3f",dihedralsXt[i]*180/M_PI);
					fprintf(stdout,"\n");

					enatRama += rama_energy2(dihedralsXt, 3, maps, size);

					// Ct anchor contribution
					finddihedral( coordCt, 6, dihedralsXt+3 ); // +3 required by "rama_energy2()"
					dihedralsXt[5] = dihedrals[nresl*3-1]; // Load the Ct anchor's Phi (mobile)

					fprintf(stdout,"Ct dihedrals: ");
					for(int i=0; i<9; i++)
						fprintf(stdout," %8.3f",dihedralsXt[i]*180/M_PI);
					fprintf(stdout,"\n");

					enatRama += rama_energy2(dihedralsXt, 3, maps, size);
					fprintf(stderr,"RAMA_ENERGY Native: %f\n",enatRama);

					enat += rama_factor * enatRama;
				}

				// fprintf(stderr,"nresl= %d\n",nresl);

				// Compute loop energies
				double ecurrent;
				for(int n=0; n<nloops; n++) // Screen loops
				{
					// Get the dihedrals array (dihedral) for supplied coordinates (co)
					finddihedral( lcoord + nresl*9*n, nresl*3, dihedrals );

					// Dump dihedrals:
//					fprintf(stderr,"dihedrals: ");
//					for(int x=0; x<nresl*3; x++)
//						fprintf(stderr," %f",dihedrals[x]*180/M_PI);
//					fprintf(stderr,"\n");

					// Ramachandran energy derived from neighbor-dependent PDFs from Dunbrack's paper.
//					energies[n] = rama_energy2(dihedrals, nresl, maps, size);
					ecurrent = rama_energy2(dihedrals, nresl, maps, size);
//					fprintf(stderr,"mobile loop energy: %f\n",ecurrent);

					// Nt anchor contribution
					finddihedral( coordNt, 6, dihedralsXt );
					dihedralsXt[6] = dihedrals[3]; // Load the Nt anchor's Psi (mobile)

//					fprintf(stderr,"Nt dihedrals: ");
//					for(int i=0; i<9; i++)
//						fprintf(stderr," %8.3f",dihedralsXt[i]*180/M_PI);
//					fprintf(stderr,"\n");

					ecurrent += rama_energy2(dihedralsXt, 3, maps, size);
//					fprintf(stderr,"Nt anchor energy: %f\n",rama_energy2(dihedralsXt, 3, maps, size));

					// Ct anchor contribution
					finddihedral( coordCt, 6, dihedralsXt+3 );
					dihedralsXt[5] = dihedrals[nresl*3-1]; // Load the Ct anchor's Phi (mobile)

//					fprintf(stderr,"Ct dihedrals: ");
//					for(int i=0; i<9; i++)
//						fprintf(stderr," %8.3f",dihedralsXt[i]*180/M_PI);
//					fprintf(stderr,"\n");

					ecurrent += rama_energy2(dihedralsXt, 3, maps, size);
//					fprintf(stderr,"Ct anchor energy: %f\n",rama_energy2(dihedralsXt, 3, maps, size));

					energies[n] += rama_factor * ecurrent; // Ramma factor included

//					fprintf(stderr,"Loop %4d RAMA_ENERGY: %f\n", n, ecurrent);

					if(n % 10 == 0)
						indicator("energy> Rama energy calculation",n,nloops);
				}
				free(dihedrals);
				free(dihedralsXt);
				free(maps);
				free(lcoord);
				indicator("energy> Rama energy calculation",nloops,nloops);
			}
			else // Single PDB file (per target) with one structure (or processing a List of PDBs)
			{
//				fprintf(stderr,"nres= %3d --> seq= %s\n",nres,seq);

				int *iaa = NULL; // AA index for each loop aminoacid
				iaa = seq2iaa(seq,nres);

//				fprintf(stderr,"iaa: ");
//				for(int i=0; i < nres; i++) // Screen  AAs (mobile residues only)
//					fprintf(stderr," %d",iaa[i]);
//				fprintf(stderr,"\n");

				float ***maps=NULL; // PDF-maps for each loop aminoacid
				if( !(maps = (float ***) malloc(sizeof(float **) * nres) ) )
				{
					printf("Sorry, unable to allocate maps memory!!!\n");
					exit(1);
				}

				//  Generating neighbor dependent distributions for current loop
				int ileft, iright; // left/right 1-letter code aminoacids
				for(int i=1; i < nres-1; i++) // Screen  AAs (mobile residues only)
				{
					// fprintf(stderr,"AA=%d  iaa[i]=%d iaa[i-1]=%d  iaa[i+1]=%d  size=%d\n",i,iaa[i],iaa[i-1],iaa[i+1],size);
					if(iaa[i-1]==20) // If cis-Pro to the left... (all Prolines are merged when Left or Right side, see Dunbracks paper...)
						ileft = 12; // left-side cis-Pro is set back to trans-Pro (code 12)
					else
						ileft = iaa[i-1];
					if	(iaa[i+1]==20) // If cis-Pro to the right... (all Prolines are merged when Left or Right side, see Dunbracks paper...)
						iright = 12; // right-side cis-Pro is set back to trans-Pro (code 12)
					else
						iright = iaa[i+1];

//					fprintf(stderr,"current residue: %c %3d --> c= %d  l= %d  r= %d\n",seq[i],i,iaa[i],ileft,iright);

//					show_map(allmaps[iaa[i]][ileft][iright],size,"deleteme.txt"); // Gnuplot command:  plot "mierda_01.txt" using 1:2:3 with image
					maps[i] = allmaps[iaa[i]][ileft][iright]; // assign the corresponding Rama-energy map to every residue
				}

				// N- or C-terminal aminoacids do not have valid Ramachandran (Nt residue lacks Phi and Ct lacks Psi)
				maps[0] = NULL;
				maps[nres-1] = NULL;
				free(iaa);

				// Allocate memory for current target dihedrals
				float *dihedrals = (float *) malloc( sizeof(float) * nres*3 );

				fprintf(stderr,"Get the dihedrals array (dihedral) for supplied coordinates (co)\n");
				// Get the dihedrals array (dihedral) for supplied coordinates (co)
				finddihedral( coord, nres*3, dihedrals );

//				fprintf(stderr,"dihedrals: ");
//				for(int x=0; x<nres*3; x++)
//					fprintf(stderr," %f",dihedrals[x]);
//				fprintf(stderr,"\n");

				// Ramachandran energy derived from neighbor-dependent PDFs from Dunbrack's paper.
//				energy = rama_energy2(dihedrals, nres, maps, size);
				energy += rama_factor * rama_energy2(dihedrals, nres, maps, size);

				free(dihedrals);
				free(maps);
			}
		}

		// "allmaps" should be freed some day...
//		if(rama_model != 0)
//			for(int i=1; i < nresl-1; i++) // Screen LOOP AAs (mobile residues only)
//			{
//				for(int j=0; j<size; j++)
//					free(maps[i][j]);
//				free(maps[i]);
//			}


		if(outfile_switch)
		{
			// sprintf(dummy,"%s_%s",pdbnames[n],outfile);
			fout = fopen(outfile,"w");
			if(fout == NULL) // problem while opening file...
			{
				fprintf(stdout,"%s> Output file %s could not be opened...\n",prog,outfile);
				exit(1);
			} // File opening check passed...
			// fprintf(fout,"# index  Energy\n");
		}

		//
		// OUTPUT
		//
		if(!loops_switch && !many_ligands)
		{
			if(npdbs > 1) // Multiple run mode
			{
				fprintf(stdout,"%s> %s energy: %f\n",prog,pdbnames[n],energy);
				all_energies[n] = energy; // Store current energy for single-run mode (computing energies from a long complete structures list)
			}
			else // Single run
			{
				if(debug)
					fprintf(stdout,"Energy: ");
				fprintf(stdout,"%s> Energy: %f\n",prog,energy);
				if(debug)
					fprintf(stdout,"\n");
			}
		}
		else // Loops and others...
		{
			if(outfile_switch)
			{
				if(rmsd_switch) // RMSD Nat
					if(rmsdfile_switch) // RMSD Ref
					{
						fprintf(fout, "#%-4s %12s %12s %12s\n","loop", "Energy", "RMSD_Nat", "RMSD_Ref");
						for(int n=0; n<nloops; n++) // Screen loops
							fprintf(fout,"%5d %12.6f %12.6f %12.6f\n",n+1,energies[n],rmsds[n],rmsds2[n]);
					}
					else
					{
						fprintf(fout, "#%-4s %12s %12s\n","loop", "Energy", "RMSD_Nat");
						for(int n=0; n<nloops; n++) // Screen loops
							fprintf(fout,"%5d %12.6f %12.6f\n",n+1,energies[n],rmsds[n]);
					}
				else
					if(rmsdfile_switch) // RMSD Ref
					{
						fprintf(fout, "#%-4s %12s %12s\n","loop", "Energy", "RMSD_Ref");
						for(int n=0; n<nloops; n++) // Screen loops
							fprintf(fout,"%5d %12.6f %12.6f\n",n+1,energies[n],rmsds2[n]);
					}
					else
					{
						fprintf(fout, "#%-4s %12s\n","loop", "Energy");
						for(int n=0; n<nloops; n++) // Screen loops
							fprintf(fout,"%5d %12.6f\n",n+1,energies[n]);
					}
			}

			if(has_native && outfile_switch)
			{
				if(rmsd_switch) // RMSD Nat
					if(rmsdfile_switch) // RMSD Ref
						fprintf(fout,"# native %5d %12.6f %12.6f %12.6f\n",0,enat,0.0,0.0);
					else
						fprintf(fout,"# native %5d %12.6f %12.6f\n",0,enat,0.0);
				else
					if(rmsdfile_switch) // RMSD Ref
						fprintf(fout,"# native %5d %12.6f %12.6f\n",0,enat,0.0);
					else
						fprintf(fout,"# native %5d %12.6f\n",0,enat);

//				if(rmsd_switch) // Computing RMSD
//					fprintf(fout,"%5d %12.6f %12.6f\n",0,enat,0.0);
//				else // Without computing RMSD (faster)
//					fprintf(fout,"%5d %12.6f\n",0,enat);
			}

			if(outfile_switch)
				fprintf(stdout,"Energies were dumped into: %s\n",outfile);

			// Swapping RMSDs for statistics
			if(rmsdfile_switch)
			{
				float *kk;
				kk = rmsds;
				rmsds = rmsds2;
				rmsds2 = kk;
			}

			// Generating GNU-Plot output
			FILE *f_plot; // Plots file
			if(gplot_switch)
			{
				// sprintf(dummy,"%s.gp",outfile);
				if( !(f_plot = fopen(gpfile,"w") ) )
				{
					printf( "%s> Error, I can't write into %s text file! Forcing exit!\n", prog, gpfile );
					exit(1);
				}

				// Ploting Energies vs. RMSD
				sprintf(dummy,"%s_%s",pdbnames[n],outbase); // Setting single output results file name
				fprintf(f_plot,"set title \"%s\"\n unset key\n ",dummy);
				fprintf(f_plot,"set terminal jpeg size 640,480\n");
				fprintf(f_plot,"set xrange [0:*]\n");
				fprintf(f_plot,"set xlabel \"RMSD [A]\"\n");
				fprintf(f_plot,"set ylabel \"Energy [a.u.]\"\n");
				fprintf(f_plot,"set output \"%s.jpg\"\n",dummy);
				fprintf(f_plot,"plot \"-\" u 1:2 w p ps 0.2, \"-\" u 1:2 w p ps 1 lt 3\n");

				for(int n=0; n<nloops; n++) // Screen loops
					fprintf(f_plot,"%10.4f %10.4f\n", rmsds[n], energies[n]);

				fprintf(f_plot,"e\n0.0 %10.4f\ne\n", enat);

				fclose(f_plot);
				fprintf(stdout,"Gnuplot energies dumped into: %s\n",gpfile);
			}

// exit(0);

			// ------------------------------------------------------------------------------------
			// Statistics computation (usually require RMSD)
			// ------------------------------------------------------------------------------------
			if(rmsd_switch)
			{
				// Sort array of energies
				int *sorte = (int *) malloc(sizeof(int) * nloops); // Allocate memory for the SORTed indices array for Energies (sorte)
				for(int u = 0; u < nloops; u++) // Initialize sorted array of energies
					sorte[u] = u; // Load with the sequential indices
				quicksort(energies, sorte, 0, nloops-1); // Quick sort algorithm
				if(debug)
					for(int i=0; i<nloops; i++)
						fprintf(stdout,"%4d %12.6f %12.6f\n", i+1, energies[sorte[i]], rmsds[sorte[i]]);

				// Energy cutoff to get the best energies
				float ecut = energies[ sorte[0] ];

				if(outfile_switch)
					fprintf(fout,"#\n# Lowest energy loop: %4d  Energy= %f  RMSD= %f\n", sorte[0]+1, ecut, rmsds[sorte[0]]);
				all_bestrmsds[n] = rmsds[sorte[0]];

				// Sort array of RMSDs
				int *sortr = (int *) malloc(sizeof(int) * nloops); // Allocate memory for the SORTed indices array for RMSDs (sortr)
				for(int u = 0; u < nloops; u++) // Initialize sorted array of RMSDs
					sortr[u] = u; // Load with the sequential indices
				quicksort(rmsds, sortr, 0, nloops-1); // Quick sort algorithm
				if(debug)
					for(int i=0; i<nloops; i++)
						fprintf(stdout,"%4d %12.6f %12.6f\n", i+1, energies[sortr[i]], rmsds[sortr[i]]);

				// Rank of best loop
				int rank = -1;
				for(int i=0; i<nloops; i++)
					if(sortr[0] == sorte[i])
					{
						rank = i+1;
						break;
					}
				if(outfile_switch)
					fprintf(fout,"# Best loop rank: %4d\n", rank);
				all_bestranks[n] = (float) rank;

				// Native Z-score
				all_zscores[n] = zscore(energies, enat, nloops);
				if(outfile_switch)
					fprintf(fout,"# Z-score (native): %f\n#\n", all_zscores[n]); // Compute the z-score of "n" data elements ("data") wrt some reference value ("ref")

				// Pearson's correlation coefficient ("r" factor)
				all_corrs[n] = corrPearson(rmsds, energies, nloops);
				if(outfile_switch)
					fprintf(fout,"# Correlation (rmsd vs. energy): %f\n#\n", all_corrs[n]);

				// Enrichment score
				if(outfile_switch)
				{
					int nbetter; // Number of cases better than current threshold
					int nbest; // Number of best cases
					float *work = (float *) malloc( sizeof(float) * nloops); // Allocate the maximum workspace required
					float irmsd; // Current RMSD threshold
					float maxrmsd; // Maximum RMSD
					maxrmsd = getmax(rmsds,nloops);

					fprintf(fout,"# Enrichment score [%%]:\n# --------------------\n# %%best"); // Dump header
					for(irmsd = drmsd; irmsd < maxrmsd; irmsd += drmsd)
						fprintf(fout," <%4.1f",irmsd);
					fprintf(fout,"\n#");

					for(int nbin = 0; nbin < nebins; nbin++)
					{
						nbest = (int) ((float)nloops*ebins[nbin]);
						fprintf(fout," %5.1f", 100*ebins[nbin]);

						// Select the rmsds for the "nbest" (lowest energy) loops
						for(int j=0; j<nbest; j++)
							work[j] = rmsds[ sorte[j] ];

						nbetter = 0;
						for(irmsd = drmsd; irmsd < maxrmsd; irmsd += drmsd)
						{
							nbetter += numinrange(work, nbest, irmsd-drmsd, irmsd);
							// fprintf(stdout,"In the best %d loops, %d are better than %3.1f A RMSD: enrichment= %3.0f\n", nbest, nbetter, irmsd, 100*(float)nbetter/nbest);
							fprintf(fout," %5.1f", 100*(float)nbetter/nbest);
						}
						fprintf(fout,"\n#");
					}

					fprintf(fout,"\n# HIST:");
					for(irmsd = drmsd; irmsd < maxrmsd; irmsd += drmsd)
					{
						fprintf(fout," %5.1f", 100*(float)numinrange(work, nbest, irmsd-drmsd, irmsd)/nloops); // Work must contain the 100%
					}
					fprintf(fout,"\n#");

					free(work);
				}

				// Compute the N-best RMSD (from fraction of the total number of loops, "nloops")
				for(int x = nbmax-1; x >= 0 ; x--) // screen backwards to account for "not enough models per case"
				{
					float bestrmsd = 999999; // Some high value
					for(int b=0; b < (int)roundf(nloops*fb[x]); b++)
						if(	rmsds[ sorte[b] ] < bestrmsd ) // If better...
							bestrmsd = rmsds[ sorte[b] ]; // ...save it!
					if(bestrmsd==999999 && x < nbmax-1) // If not enough models per case...
						all_nbestrmsds[x][n] = all_nbestrmsds[x+1][n]; // ...and dump it!
					else // Normal behavior
						all_nbestrmsds[x][n] = bestrmsd; // ...and dump it!
				}

				all_rmsds[n] = rmsds; // Store current case RMSDs into "all_rmsds" array
				all_sorte[n] = sorte; // Store current case Energy Sorted indices into "all_sorte" array
				all_nloops[n] = nloops; // Store current numer of loops into "all_nloops" array

				// Free stuff
				free(sortr); // not needed anymore...
				// free(sorte);
				free(energies);
				// free(rmsds);

			}

			// Close output file
			if(outfile_switch)
				fclose(fout);
		}

		// Not needed anymore?
		if(emodel == 5 && map->frame_model == 10 && mon_switch ) // Only required for KORP
		{
			free(resnums); // Main-PDB residue numbers (loop-less, i.e. environment)
			free(reschainids); // Main-PDB chain-ids (loop-less, i.e. environment)
			free(resnumsnat); // Main-PDB Native residue numbers
			free(reschainidsnat); // Main-PDB Native chain-ids
			free(resnumsloop); // Loop residue numbers (only mobile)
			free(reschainidsloop); // Loop chain-ids (only mobile)
		}
	}

	// RMSD profile
	if(npdbs > 1) // Multiple run mode
	{
		sprintf(dummy,"%s_score.txt",outbase); // Setting single output results file name

		FILE *fscore;

		if(loops_switch && rmsd_switch) // loops and RMSD requested
		{
			fscore = fopen(dummy,"w");
			if(fscore == NULL) // problem while opening file...
			{
				fprintf(stdout,"%s> Output file %s could not be opened...\n",prog,dummy);
				exit(1);
			} // File opening check passed...

			float dprof = 0.2; // Delta RMSD for RMSD profile generation
			fprintf(fscore,"%5s|[%%]","#RMSD");
			for(int x = 0; x < nbmax; x++)
				fprintf(fscore," %4.1f",100.0*fb[x]);
			fprintf(fscore,"\n");

			for(float r = 0.0; r <= 5.0; r += dprof)
			{
				fprintf(fscore,"%5.1f    ",r);
				for(int x = 0; x < nbmax; x++)
				{
					int nbetter = 0; // Number of cases where at least one of the "nb" lowest energy loops has RMSD below "r" threshold
					for(int n=0; n < npdbs; n++)
					{
						for(int b=0; b < (int)roundf(all_nloops[n] * fb[x]); b++)
							if(	all_rmsds[n][ all_sorte[n][b] ] < r )
							{
								nbetter++;
								break;
							}
					}
					fprintf(fscore," %4.0f", 100*(float)nbetter/(float)npdbs );
				}
				fprintf(fscore,"\n");
			}

			for(int n=0; n < npdbs; n++)
			{
				// Not needed anymore
				free(all_rmsds[n]);
				free(all_sorte[n]);
			}
			free(all_rmsds);
			free(all_sorte);

			float avg_zscores = average(all_zscores,npdbs); // All Z-scores array (for statistics computation)
			float avg_bestrmsds = average(all_bestrmsds,npdbs); // All best-RMSDs array (for statistics computation)
			float avg_bestranks = average(all_bestranks,npdbs); // All best-Ranks array (for statistics computation)
			float avg_corrs = average(all_corrs,npdbs); // All cross-correlations array (for statistics computation)
			fprintf(fscore,"#\n#---------------------------------------------------------- %%Best_RMSD [A]:\n");
			fprintf(fscore,"%-20s %7s %9s %9s %9s ", "#PDB-ID", "Z-score", "LowE-RMSD", "LowE-Rank", "CrossCorr");
			for(int x = 0; x < nbmax; x++)
				fprintf(fscore," %5.1f%%", fb[x]*100);
			fprintf(fscore,"\n");
			fprintf(fscore,"#----------------------------------------------------------");
			for(int x = 0; x < nbmax; x++)
				fprintf(fscore,"%s", "-------");
			fprintf(fscore,"\n");
			for(int n=0; n<npdbs; n++)
			{
				fprintf(fscore,"%-20s %7.3f %9.3f %9.0f %9.3f ", pdbnames[n], all_zscores[n], all_bestrmsds[n],
						all_bestranks[n], all_corrs[n]);

				for(int x = 0; x < nbmax; x++)
					fprintf(fscore," %6.3f", all_nbestrmsds[x][n]); // ...and dump it!
				fprintf(fscore,"\n");
			}
			fprintf(fscore,"#----------------------------------------------------------");
			for(int x = 0; x < nbmax; x++)
				fprintf(fscore,"%s", "-------");
			fprintf(fscore,"\n");

			// Computing general averages
			fprintf(fscore,"%-20s %7.3f %9.3f %9.3f %9.3f ","#AVERAGE:", avg_zscores, avg_bestrmsds, avg_bestranks, avg_corrs );
			for(int x = 0; x < nbmax; x++)
				fprintf(fscore," %6.3f", average(all_nbestrmsds[x],npdbs)); // ...and dump it!
			fprintf(fscore,"\n");

			// Computing general sigmas (standard deviations)
			fprintf(fscore,"%-20s %7.3f %9.3f %9.3f %9.3f ","#SIGMA:", sigma(all_zscores,npdbs,avg_zscores),
					sigma(all_bestrmsds,npdbs,avg_bestrmsds), sigma(all_bestranks,npdbs,avg_bestranks),
					sigma(all_corrs,npdbs,avg_corrs) );
			for(int x = 0; x < nbmax; x++)
				fprintf(fscore," %6.3f", sigma(all_nbestrmsds[x],npdbs)); // ...and dump it!
			fprintf(fscore,"\n");

			// Computing general medians
			fprintf(fscore,"%-20s %7.3f %9.3f %9.3f %9.3f ","#MEDIAN:", median(all_zscores,npdbs),
					median(all_bestrmsds,npdbs), median(all_bestranks,npdbs),
					median(all_corrs,npdbs) );
			for(int x = 0; x < nbmax; x++)
				fprintf(fscore," %6.3f", median(all_nbestrmsds[x],npdbs)); // ...and dump it!
			fprintf(fscore,"\n");

			free(all_zscores);
			free(all_bestrmsds);
			free(all_bestranks);
			free(all_corrs);
			for(int x = 0; x < nbmax; x++)
				free(all_nbestrmsds[x]);
			free(all_nbestrmsds);

			fclose(fscore);
			fprintf(stdout,"%s> Output summary of scores dumped into %s\n",prog,dummy);
		}

		if(!loops_switch) // protein modeling
		{
			fscore = fopen(dummy,"w");
			if(fscore == NULL) // problem while opening file...
			{
				fprintf(stdout,"%s> Output file %s could not be opened...\n",prog,dummy);
				exit(1);
			} // File opening check passed...

			fprintf(fscore,"%-20s %12s\n", "#PDB-name", "KORP_energy");
			// fprintf(fscore,"#--------------------------------\n");
			for(int n=0; n<npdbs; n++)
				fprintf(fscore,"%-20s %12.6f\n", pdbnames[n], all_energies[n]);

			fclose(fscore);
			fprintf(stdout,"%s> Output summary of energies dumped into %s\n",prog,dummy);
		}


//		free(all_zscores);
//		free(all_bestrmsds);
//		free(all_bestranks);
//		free(all_corrs);
//		for(int x = 0; x < nbmax; x++)
//			free(all_nbestrmsds[x]);
//		free(all_nbestrmsds);
	}

}
// #########################################################################
// END OF MAIN
// #########################################################################



// #########################################################################
// FUNCTIONS ZONE
// #########################################################################

