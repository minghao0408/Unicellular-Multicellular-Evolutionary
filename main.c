// ************* START SECTION - Credits, compiling and executing instructions ************* //


// To compile:
// gcc main.c -o main -lm -Wall

// To run (Example arguments):
// ./main 200 200 200000 1000 0 100 1 5
// Args: rows, cols, duration, savestep, start_J, start_pop%, mutation_on, repr_thresh

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <string.h>

// Forward declarations of helper functions
unsigned int rand_interval(unsigned int min, unsigned int max, long long int *randcount);
void shuffle( int **array,  int n, long long int *randcount);
void thoroidal(int *row_th, int *col_th, int fieldrows, int fieldcols);
double randn (double mu, double sigma, long long int *randcount);

// Global File Pointers (Only the essential files for reproduction are retained)
FILE *Output_File;      // Records the spatial distribution of cells (0=empty, 1=cell)
FILE *Adhspec_File;     // Records the adhesion values (J) of cells
FILE *Time_File;        // Logs the start and end time of execution

int main(int argc, char **argv)
{
    fprintf(stderr,"Start Simulation...\n");
    
    // Scaling factor for energy units. 
    // Note: In the code, rho = 1 and xi = 0.01, implying a scaling factor of 100.
    float rho_xi_scalingfactor=1; 
    float eul_dt=1; // Euler integration time step
    
    // --- 1. Parse Command Line Arguments ---
    int fieldrows =atoi(argv[1]);
    int fieldcols =atoi(argv[2]);
    int duration=atoi(argv[3]);     // Total simulation steps (e.g., 200000)
    int savestep= atoi(argv[4]);    // Save data every N steps (e.g., 1000)
    float start_ww_int=atof(argv[5]); // Initial adhesion value (usually 0)
    int pop1percentage=atoi(argv[6]); // Initial population density (%)
    int adhmat_mut_spec1=atoi(argv[7]); // Allow mutation (1=yes, 0=no)
    float reproductionthresh=atof(argv[8])/rho_xi_scalingfactor; // Division energy threshold
    
    // --- 2. Simulation Variables Definition ---
    long long int randcount=0;
    double v, vmut;
    float mutatingsigma=0.1; // Mutation standard deviation (Note: Code uses 0.1, Paper says 0.005)
    float mutationprob=0.005; // Probability of mutation per division
    int startsem;
    int actualPatchType=-1000;
    
    // Vectors for shuffling positions (to ensure random update order)
    int posvectlength = (fieldrows-1)*(fieldcols-1); 
    int neighbvectlength = 9;
    int neighbrows = 3;
    int neighbcols = 3;
    
    // Temporary variables for neighbors
    int neighbPatchType=-1000;
    int na_patchtype=-1000;
    int nc_patchtype=-1000;
    
    // Counters
    int countersavestep=0; 
    int numberofcells=0;
    long long int integralnumberofcells=0;
    
    double TT, dt, tRunning;
    tRunning = 0.;
    dt=1;
    
    // --- 3. Biological Rules ---
    // Death parameters
    int randomdeathallowed=1; 
    float dierandthresh=0.1; // Random death probability per step
    int agedeathallowed=0;   // Aging is DISABLED in this reproduction version
    int maximumAge=50;
    int metabolismallowed=1; 
    int deathallowed=metabolismallowed; 
    int energydeathallowed = 1; // Death by starvation allowed
    int birthallowed=1; 
    
    // Diffusion and Physics Parameters
    int difftimemult=1; // How many diffusion steps per Monte Carlo step
    float dierand;
    float T = 1; // Temperature for Metropolis algorithm (noise)
    float fooddiffpercent=0.05/difftimemult; // Diffusion coefficient D_R
    float fooddiffpercentcell=0.05/difftimemult; 
    
    // Energy parameters (Scaled)
    float minstartenergy=1/rho_xi_scalingfactor; 
    float randstartenergy=9/rho_xi_scalingfactor; 
    int absorbedenergycellcapping=1; 
    float maxcellabsorbedenergy=reproductionthresh/rho_xi_scalingfactor; 
    float death_lacken_thresh=1/rho_xi_scalingfactor; // Starvation threshold (B_min)
    float transienttime=20; // Time before inoculation
    float Su_val;
    
    // Metabolism Parameters
    float subst_t0=17.0; // Initial resource level
    float phi_par=0.085; // Resource input rate
    float eta_N_par=0.005; // Resource decay rate
    float xi_par=0.01;   // Resource consumption rate (Note: 0.01 in code vs 0.001 in paper)
    float eta_E_par=0.1; // Metabolic cost
    float rho_par=1/rho_xi_scalingfactor; // Energy conversion efficiency (1.0)
    
    // Temporary calculation variables
    double potential1, potential2;
    float H, tempEn, tempSu;
    float tempAdh;
    int tempAgeMap;
    float tempBiomassMap;
    float rn;
    float Ds, Dscell;
    
    // --- 4. File Initialization ---
    Output_File = fopen("outputfile.txt", "w");
    Adhspec_File = fopen("adhspecfile.txt", "w");
    
    // Random Seed Initialization
    startsem=1551136312; // Fixed seed used in the original paper for reproducibility
    srand(startsem); 
    
    // Loop iterators
    int neighbcontainszero, selected, notme;
    int i,j, ii, jj, irow, icol, nrow, ncol, i_radius, j_radius, napatchrow, napatchcol, ncpatchrow, ncpatchcol, fi, fj, iii, jjj, poscount, pos, neighbposcount, neighbpos, difftimecount;
    int tempXdiffpos, tempYdiffpos;
    
    int foodpatchradius;
    int inoculationradius;
    
    // --- 5. Simulation Loop Control ---
    // Modify these values to change the Resource Size (phi)
    int startingvarfoodpatch=45; // STARTING Radius (phi/2)
    int deltavarfoodpatch=4;     // Step size (if scanning multiple sizes)
    int maxfoodpatch=45;         // MAX Radius (Set equal to startingvarfoodpatch for single run)
    int maxstatisticsiteration=1;// Replicates per size
    
    // --- 6. Memory Allocation ---
    // Allocating 2D arrays for the lattice grids
    float** En=(float**)calloc(fieldrows, sizeof(float *)); // Internal Energy (Biomass)
    float** Su=(float**)calloc(fieldrows, sizeof(float *)); // Substrate (Resource)
    float** newSu=(float**)calloc(fieldrows, sizeof(float *)); // Temp Substrate for diffusion
    int** patchesMap=(int**)calloc(fieldrows, sizeof(int *)); // Cell locations (0 or 1)
    float** adhesionMap=(float**)calloc(fieldrows, sizeof(float *)); // Adhesion values (J)
    int** ageMap=(int**)calloc(fieldrows, sizeof(int *)); // Cell Age
    float** biomassMap=(float**)calloc(fieldrows, sizeof(float *)); // Accumulated Biomass
    
    // Allocating helper vectors
    int** neighbMatrix=(int**)calloc(neighbrows, sizeof(int *)); 
    int** posrowscols=(int**)calloc(posvectlength, sizeof(int *)); 
    int** neighbrowscols=(int**)calloc(neighbvectlength, sizeof(int *));

    for(i=0;i<fieldrows;i++){
        En[i]=calloc(fieldcols, sizeof(float));
        Su[i]=calloc(fieldcols, sizeof(float));
        newSu[i]=calloc(fieldcols, sizeof(float));
        patchesMap[i]=calloc(fieldcols, sizeof(int));
        adhesionMap[i]=calloc(fieldcols, sizeof(float));
        ageMap[i]=calloc(fieldcols, sizeof(int));
        biomassMap[i]=calloc(fieldcols, sizeof(float));
    }
    
    for(ii=0;ii<neighbrows;ii++){
        neighbMatrix[ii]=calloc(neighbcols, sizeof(float));
    }
    
    for(i=0;i<neighbvectlength;i++){
        neighbrowscols[i]=(int*)calloc(2, sizeof(int));
    }
    
    for(i=0;i<posvectlength;i++){
        posrowscols[i]=(int*)calloc(2, sizeof(int));
    }
    
    // Log Start Time
    time_t timer;
    char buffer[26];
    struct tm* tm_info;
    timer = time(NULL);
    tm_info = localtime(&timer);
    strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);
    Time_File = fopen("timefile.txt", "w");
    fprintf(Time_File,"Start: %s \n", buffer);

    // ************ MAIN SIMULATION LOOPS ************
    // Outer loops control the experimental parameters (Resource Size)
    
    for(int varfoodpatch=startingvarfoodpatch; varfoodpatch<maxfoodpatch+1; varfoodpatch=varfoodpatch+deltavarfoodpatch){
        fprintf(stderr,"Running simulation for Foodpatch Radius: %d\n", varfoodpatch);
        
        for(int statisticsiteration=0; statisticsiteration<maxstatisticsiteration; statisticsiteration++){
            
            foodpatchradius=varfoodpatch; 
            inoculationradius=varfoodpatch; 
            
            countersavestep=0;
            numberofcells=0; 
            integralnumberofcells=0; 
            
            // A. Reset Lattice
            for(i=0;i<fieldrows;i++){ 
                for (j=0;j<fieldcols;j++){
                    patchesMap[i][j]=0;
                    adhesionMap[i][j]=0;
                    ageMap[i][j]=-1000; 
                    biomassMap[i][j]=0;
                }
            }
            
            // B. Initialize Resource Field (Circle shape)
            for(i=1;i<fieldrows;i++){
                for (j=1;j<fieldcols;j++){
                    i_radius=i-fieldrows/2;
                    j_radius=j-fieldcols/2;
                    // Place initial resource only inside the radius
                    if(sqrt(pow(i_radius,2)+pow(j_radius,2))<foodpatchradius){
                        Su[i][j]=subst_t0; 
                    }
                }
            }

            // C. Initialize Position Vectors (excluding borders for safety)
            poscount=0;
            for(i=1;i<fieldrows-1;i++){ 
                for(j=1;j<fieldcols-1;j++){ 
                    posrowscols[poscount][0]=i; 
                    posrowscols[poscount][1]=j; 
                    poscount++;
                }
            }
            
            neighbposcount=0;
            for(i=0;i<3;i++){
                for(j=0;j<3;j++){
                    neighbrowscols[neighbposcount][0]=i; 
                    neighbrowscols[neighbposcount][1]=j; 
                    neighbposcount++;
                }
            }
            
            // ------------------ TIME LOOP (Monte Carlo Steps) ------------------
            for (TT=0.0; TT<duration+1; TT+=dt){  
                
                // 1. Inoculation (Seed initial cells at transient time)
                if (TT==transienttime){
                    for(i=1;i<fieldrows;i++){
                        for (j=1;j<fieldcols;j++){
                            i_radius=i-fieldrows/2;
                            j_radius=j-fieldcols/2;
                            // Only inoculate inside the resource patch
                            if(sqrt(pow(i_radius,2)+pow(j_radius,2))<inoculationradius){
                                v = ((double)rand() / (double)RAND_MAX) * 100.0 ; 
                                randcount++;
                                
                                if( v <= pop1percentage){ // Probability check
                                    numberofcells++;
                                    patchesMap[i][j]=1; // Create cell
                                    adhesionMap[i][j]=start_ww_int; // Initial Adhesion (J=0)
                                    ageMap[i][j]=0; 
                                    biomassMap[i][j]=0;
                                    v = (double)rand() / (double)RAND_MAX ;
                                    randcount++;
                                    En[i][j]=minstartenergy+randstartenergy*v; // Random initial energy
                                }
                            }
                        }
                    }
                }
                
                // --- MODULE 1: PHYSICS (Metropolis Algorithm for Cell Movement) ---
                // We shuffle the update order to avoid bias
                shuffle(posrowscols,posvectlength, &randcount);
                
                for (pos=0;pos<posvectlength;pos++){ 
                    i=posrowscols[pos][0]; 
                    j=posrowscols[pos][1]; 
                    actualPatchType=patchesMap[i][j];
                    
                    // Only process if there is a cell at this location
                    if (actualPatchType==1){
                        // Select a random neighbor that is NOT self
                        notme=0;
                        while (notme==0){
                            nrow=rand_interval(0,3, &randcount); 
                            ncol=rand_interval(0,3, &randcount); 
                            if ( nrow * ncol !=1) { notme=1; }
                        }
                        nrow=i+nrow-1; 
                        ncol=j+ncol-1;
                        // Handle periodic boundary conditions (Thoroidal world)
                        thoroidal(&nrow,&ncol,fieldrows,fieldcols);
                        
                        neighbPatchType=patchesMap[nrow][ncol]; 
                        
                        // --- Calculate Hamiltonian Energy AFTER Swap (potential1) ---
                        potential1=0.0;
                        // Calculate energy of neighbor cell (moved to current pos)
                        for (irow = -1; irow < 2; irow++){ 
                            for (icol = -1; icol < 2; icol++){
                                napatchrow=i+irow;
                                napatchcol=j+icol;
                                thoroidal(&napatchrow,&napatchcol,fieldrows,fieldcols); 
                                if (!((napatchrow==i)&&(napatchcol==j))) {  
                                    na_patchtype=patchesMap[napatchrow][napatchcol]; 
                                    // Interaction Energy = (J_1 + J_2) / 2
                                    if(neighbPatchType==1 && na_patchtype==1){ 
                                        potential1=potential1+((adhesionMap[nrow][ncol]+adhesionMap[napatchrow][napatchcol])/2);
                                    } 
                                }
                            }
                        }
                        // Calculate energy of current cell (moved to neighbor pos)
                        for(irow = -1; irow < 2; irow++){
                            for (icol = -1; icol < 2; icol++){
                                ncpatchrow=nrow+irow;
                                ncpatchcol=ncol+icol;
                                thoroidal(&ncpatchrow,&ncpatchcol,fieldrows,fieldcols);
                                if (!((ncpatchrow==nrow)&&(ncpatchcol==ncol))) {  
                                    nc_patchtype=patchesMap[ncpatchrow][ncpatchcol]; 
                                    if(actualPatchType==1 && nc_patchtype==1){ 
                                        potential1=potential1+((adhesionMap[i][j]+adhesionMap[ncpatchrow][ncpatchcol])/2);
                                    } 
                                }
                            }
                        }
                        
                        // --- Calculate Hamiltonian Energy BEFORE Swap (potential2) ---
                        potential2=0.0;
                        // Energy of current cell at current pos
                        for (irow = -1; irow < 2; irow++){
                            for (icol = -1; icol < 2; icol++){
                                napatchrow=i+irow;
                                napatchcol=j+icol;
                                thoroidal(&napatchrow,&napatchcol,fieldrows,fieldcols);
                                if (!((napatchrow==i)&&(napatchcol==j))) {  
                                    na_patchtype=patchesMap[napatchrow][napatchcol]; 
                                    if(actualPatchType==1 && na_patchtype==1){ 
                                        potential2=potential2+((adhesionMap[i][j]+adhesionMap[napatchrow][napatchcol])/2);
                                    } 
                                }
                            }
                        }
                        // Energy of neighbor cell at neighbor pos
                        for (irow = -1; irow < 2; irow++){
                            for (icol = -1; icol < 2; icol++){
                                ncpatchrow=nrow+irow;
                                ncpatchcol=ncol+icol;
                                thoroidal(&ncpatchrow,&ncpatchcol,fieldrows,fieldcols);
                                if (!((ncpatchrow==nrow)&&(ncpatchcol==ncol))) {  
                                    nc_patchtype=patchesMap[ncpatchrow][ncpatchcol]; 
                                    if(neighbPatchType!=0 && nc_patchtype!=0){ 
                                        potential2=potential2+((adhesionMap[nrow][ncol]+adhesionMap[ncpatchrow][ncpatchcol])/2);
                                    } 
                                }
                            }
                        }
                        
                        // Calculate Delta H
                        H = potential1-potential2;
                        rn = (double)rand() / (double)RAND_MAX ;
                        randcount++;
                        
                        // Metropolis Acceptance Rule
                        // Accept if Energy decreases OR with probability exp(-DH/T)
                        if (H<0 || (H>=0 && rn <= (exp(-H/T)))){
                            // Execute Swap: Swap Type, Energy, Substrate, Adhesion, Age, Biomass
                            patchesMap[i][j]=neighbPatchType;
                            tempEn=En[i][j]; 
                            tempSu=Su[i][j]; 
                            tempAdh=adhesionMap[i][j];
                            adhesionMap[i][j]=adhesionMap[nrow][ncol];
                            adhesionMap[nrow][ncol]=tempAdh;
                            En[i][j]=En[nrow][ncol];
                            patchesMap[nrow][ncol]=actualPatchType;
                            En[nrow][ncol]=tempEn;
                            
                            tempAgeMap=ageMap[i][j];
                            tempBiomassMap=biomassMap[i][j];
                            ageMap[i][j]=ageMap[nrow][ncol];
                            ageMap[nrow][ncol]=tempAgeMap;
                            biomassMap[i][j]=biomassMap[nrow][ncol];
                            biomassMap[nrow][ncol]=tempBiomassMap;
                        }
                    } 
                } 
                
                // --- MODULE 2: METABOLISM (Growth, Death, Division) ---
                shuffle(posrowscols,posvectlength, &randcount);
                
                for (pos=0;pos<posvectlength;pos++){ 
                    i=posrowscols[pos][0]; 
                    j=posrowscols[pos][1]; 
                    int dead=0;
                    actualPatchType=patchesMap[i][j];
                    
                    // Define Resource Field Input
                    // Resources are only input inside the foodpatch radius
                    i_radius=i-fieldrows/2;
                    j_radius=j-fieldcols/2;
                    if(sqrt(pow(i_radius,2)+pow(j_radius,2))<foodpatchradius){
                        Su_val=phi_par;
                    }else{
                        Su_val=0.0;
                    }
                    
                    if (actualPatchType==0){
                        // Medium: Resource input and natural decay
                        En[i][j]=0;
                        newSu[i][j]=newSu[i][j] + eul_dt * (Su_val - (eta_N_par * Su[i][j])); 
                    }else if (actualPatchType==1){
                        // Cell: Metabolism
                        ageMap[i][j]= ageMap[i][j]+1;
                        
                        // Growth Equation: Gain Energy - Metabolic Cost
                        En[i][j] = En[i][j] + eul_dt * (- eta_E_par * En[i][j] + ( rho_par * Su[i][j]) ); 
                        
                        // Resource Consumption Equation
                        // Note: Cells consume resources proportional to xi_par
                        if (absorbedenergycellcapping==1 && En[i][j] >= maxcellabsorbedenergy){
                            En[i][j]=maxcellabsorbedenergy;
                            newSu[i][j] = newSu[i][j] + eul_dt * (Su_val - (eta_N_par + xi_par ) * Su[i][j]); 
                        }else{
                            newSu[i][j] = newSu[i][j] + eul_dt * (Su_val - (eta_N_par + xi_par ) * Su[i][j]); 
                        }
                        biomassMap[i][j]=biomassMap[i][j]+En[i][j]; 
                        
                        // Death Logic
                        dierand=(double)rand() / (double)RAND_MAX ; 
                        randcount++;
                        // Die if: Energy too low OR Random chance OR Too old (disabled)
                        if (((En[i][j] < death_lacken_thresh && energydeathallowed) || (dierand<dierandthresh && randomdeathallowed) || (ageMap[i][j]>=maximumAge && agedeathallowed)) && deathallowed ){
                            
                            patchesMap[i][j]=0;
                            En[i][j]=0;
                            ageMap[i][j]=-1000;
                            biomassMap[i][j]=0;
                            adhesionMap[i][j]=0;
                            dead=1;
                            numberofcells=numberofcells-1;
                            actualPatchType=0;
                        }
                        
                        if (dead!=1){ 
                            actualPatchType=patchesMap[i][j]; 
                            
                            // Division Logic
                            if (En[i][j] >= reproductionthresh && birthallowed){
                                // Check if there is space to reproduce (Moored neighborhood)
                                neighbcontainszero=0;  
                                for (ii = 0; ii < 3; ii++) {
                                    for (jj = 0; jj < 3; jj++) {
                                        neighbMatrix[ii][jj] = patchesMap[i + ii -1][j +jj -1];
                                        if (neighbMatrix[ii][jj]==0){
                                            neighbcontainszero=1; 
                                        }
                                    }
                                }
                                if (neighbcontainszero==1){ 
                                    // Find a specific empty spot
                                    selected=0;
                                    neighbpos=0;
                                    shuffle(neighbrowscols,neighbvectlength, &randcount);
                                    while (selected==0 && neighbpos<9){
                                        nrow=neighbrowscols[neighbpos][0]; 
                                        ncol=neighbrowscols[neighbpos][1]; 
                                        
                                        if (neighbMatrix[nrow][ncol]==0){
                                            nrow=i+nrow-1; 
                                            ncol=j+ncol-1;
                                            // Handle borders
                                            if (nrow==0 || nrow==fieldrows-1 || ncol ==0  || ncol==fieldcols-1 ) {
                                                thoroidal(&nrow,&ncol,fieldrows,fieldcols);
                                                if (patchesMap[nrow][ncol]==0){
                                                    selected=1;
                                                }
                                            }
                                            else{
                                                selected=1;
                                            }
                                        }
                                        neighbpos++;
                                    }
                                    
                                    // Execute Reproduction
                                    if (selected==1){ 
                                        numberofcells++;
                                        patchesMap[nrow][ncol]=actualPatchType;
                                        
                                        // Evolutionary Step: Adhesion Mutation
                                        v = randn(0,mutatingsigma,&randcount); 
                                        vmut = (double)rand() / (double)RAND_MAX ;
                                        randcount++;
                                        
                                        if (adhmat_mut_spec1==1 && vmut<=mutationprob){
                                            // Mutation occurs: J_new = J_old + random_gauss
                                            adhesionMap[nrow][ncol]=adhesionMap[i][j]+v;
                                        }else{
                                            // Inheritance
                                            adhesionMap[nrow][ncol]=adhesionMap[i][j];
                                        }
                                        
                                        // Split Energy between mother and daughter
                                        En[nrow][ncol]=En[i][j]/2;
                                        En[i][j]=En[i][j]/2;
                                        
                                        // Reset state for new cell
                                        ageMap[i][j]=0; 
                                        ageMap[nrow][ncol]=0; 
                                        biomassMap[i][j]=0;
                                        biomassMap[nrow][ncol]=0;
                                    }
                                }
                            }
                        }
                    }
                } 
                
                integralnumberofcells = integralnumberofcells+numberofcells;
                countersavestep++;
                
                // --- MODULE 3: SAVING DATA ---
                
                // Only save data at specific intervals (savestep)
                if (TT==0 || countersavestep==savestep){ 
                    if (varfoodpatch==startingvarfoodpatch && statisticsiteration==0){
                        fprintf(stderr,"Step: %.0f, Cells: %d\n", TT, numberofcells);
                        int pi,pj;
                        for(pi=0;pi<fieldrows;pi++) {
                            for(pj=0;pj<fieldcols;pj++) {
                                // Save cell positions and adhesion values to CSV-like format
                                if(pj==fieldcols-1){
                                    fprintf(Output_File,"%d",patchesMap[pi][pj]);
                                    fprintf(Adhspec_File,"%f",adhesionMap[pi][pj]);
                                }
                                else{
                                    fprintf(Output_File,"%d,",patchesMap[pi][pj]);
                                    fprintf(Adhspec_File,"%f,",adhesionMap[pi][pj]);
                                }
                            }
                            fprintf(Output_File,"\n");
                            fprintf(Adhspec_File,"\n");
                        }
                    }
                    countersavestep=0;
                }
                
                tRunning += dt;             
                
                // --- MODULE 4: DIFFUSION ---
                
                // Diffusion Logic (Executed 'difftimemult' times per step)
                for(difftimecount=0;difftimecount<difftimemult;difftimecount++){
                    for(fi=1;fi<fieldrows-1;fi++){ 
                        for(fj=1;fj<fieldcols-1;fj++){ 
                            Ds=fooddiffpercent*Su[fi][fj];
                            Dscell=fooddiffpercentcell*Su[fi][fj];
                            for( iii=-1;iii<2;iii++){
                                for (jjj=-1;jjj<2;jjj++){
                                    
                                    // Moore Neighborhood Check (8 Neighbors)
                                    // Condition ((abs(iii)+abs(jjj))!=0) ensures we check all 8 neighbors including diagonals
                                    if ((abs(iii)+abs(jjj))!=0){ 
                                        
                                        tempXdiffpos=fi+iii;
                                        tempYdiffpos=fj+jjj;
                                        
                                        // Standard Diffusion (Center loses Ds, Neighbor gains Ds)
                                        if ((tempXdiffpos!=0) && (tempXdiffpos!=(fieldrows-1)) && (tempYdiffpos!=0) && (tempYdiffpos!=(fieldcols-1))){
                                            if((patchesMap[fi][fj]==0) && (patchesMap[tempXdiffpos][tempYdiffpos]==0)){ 
                                                newSu[fi][fj]=newSu[fi][fj]-Ds; 
                                                newSu[tempXdiffpos][tempYdiffpos]=newSu[tempXdiffpos][tempYdiffpos]+Ds; 
                                            }
                                            else{ 
                                                newSu[fi][fj]=newSu[fi][fj]-Dscell; 
                                                newSu[tempXdiffpos][tempYdiffpos]=newSu[tempXdiffpos][tempYdiffpos]+Dscell; 
                                            }
                                        }else{
                                            // Handle borders (Thoroidal)
                                            thoroidal(&tempXdiffpos,&tempYdiffpos,fieldrows,fieldcols);
                                            if((patchesMap[fi][fj]==0) && (patchesMap[tempXdiffpos][tempYdiffpos]==0)){
                                                newSu[fi][fj]=newSu[fi][fj]-Ds; 
                                                newSu[tempXdiffpos][tempYdiffpos]=newSu[tempXdiffpos][tempYdiffpos]+Ds; 
                                            }
                                            else{ 
                                                newSu[fi][fj]=newSu[fi][fj]-Dscell; 
                                                newSu[tempXdiffpos][tempYdiffpos]=newSu[tempXdiffpos][tempYdiffpos]+Dscell; 
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    // Update Substrate Matrix
                    for (fi=1;fi<fieldrows-1;fi++){
                        for (fj=1;fj<fieldcols-1;fj++){
                            Su[fi][fj]=newSu[fi][fj]+Su[fi][fj];
                            newSu[fi][fj]=0;
                        }
                    }
                } 
                
            } // End Time Loop
        } // End Stats Loop
    } // End Varfoodpatch Loop
    
    // Final Cleanup
    timer = time(NULL);
    tm_info = localtime(&timer);
    strftime(buffer, 26, "%Y-%m-%d %H:%M:%S", tm_info);
    fprintf(Time_File,"End: %s \n", buffer);
    fclose(Time_File);
    
    fclose(Output_File);
    fclose(Adhspec_File);
    
    return 0;
    
} // END Main

// ************* Helper Functions *************

// Returns a random integer between min and max (inclusive of min, exclusive of max bucket limit)
unsigned int rand_interval(unsigned int min, unsigned int max, long long int *randcount)
{
    int r;
    const unsigned int range = max - min;
    const unsigned int buckets = RAND_MAX / range;
    const unsigned int limit = buckets * range;
    
    do
    {
        r = rand();
        *randcount=*randcount+1;
    } while (r >= limit);
    
    return min + (r / buckets);
}

// Fisher-Yates shuffle algorithm for randomizing simulation order
void shuffle(int **array, int n, long long int *randcount)
{
    if (n > 1)
    {
        size_t i;
        for (i = 0; i < n - 1; i++)
        {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
            *randcount=*randcount+1;
            int t1 = array[j][0];
            int t2 = array[j][1];
            array[j][0] = array[i][0];
            array[j][1] = array[i][1];
            array[i][0] = t1;
            array[i][1] = t2;
        }
    }
}

// Generates Gaussian random numbers (Box-Muller transform)
double randn (double mu, double sigma, long long int *randcount)
{
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;
    
    if (call == 1)
    {
        call = !call;
        return (mu + sigma * (double) X2);
    }
    
    do
    {
        U1 = -1 + ((double) rand () / RAND_MAX) * 2;
        U2 = -1 + ((double) rand () / RAND_MAX) * 2;
        W = pow (U1, 2) + pow (U2, 2);
        *randcount=*randcount+2;
    }
    while (W >= 1 || W == 0);
    
    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;
    
    call = !call;
    return (mu + sigma * (double) X1);
}

// Handles periodic boundary conditions (Wrap-around coordinates)
void thoroidal(int *row_th, int *col_th, int fieldrows, int fieldcols)
{
    if (*row_th==0){
        *row_th=fieldrows-2;
    } else if (*row_th==fieldrows-1){
        *row_th=1;
    }
    if (*col_th==0){
        *col_th=fieldcols-2;
    } else if (*col_th==fieldcols-1){
        *col_th=1;
    }
}