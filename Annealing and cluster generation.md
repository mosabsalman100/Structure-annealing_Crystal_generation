# Structure-annealing_Crystal_generation

//for thsi part apply the recenter_fix  command for lammps  in the "input.voronoi_ref_fix_recenter" .
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#define A 18//1
#define B 16 //1
#define C 14 //2
#define MASS 50.939999999999998
#define NR 2


#define N_STEP 1 	// sampling number over the timing
#define D_STEP 1	// difference in neighboring timings
#define TEMP 30
#define LAMMPS_EXE "mpirun -np 4 /home/bin/lammps_voronoi"			// lammps execution command
#define INP_1_LAMMPS_REF       "/home/mosab/resources/annealing/input_step1.relaxing"// reference input file (to be read and then be modified)(here we introduced voroni method to get the number of Interstitial and Vacancy and know the TDE 
#define INP_2_LAMMPS_REF       "/home/mosab/resources/annealing/input_step2.relaxing"// reference input file (to be read and then be modified)(here we introduced voroni method to get the number of Interstitial and Vacancy and know the TDE 
#define NLINE_INP_1_LAMMPS_REF 30	
#define NLINE_INP_2_LAMMPS_REF 30	
#define N_LINE_OUTPUT_FILE	1097				// the number of lines in the output file after first step input file
#define READ_OUTPUT "/home2/mosab/2017.01.04_V_MO/2017_12_12_V_STRAIN/annealing/output.voronoi_step_01"	// restart file name

#define READ_FILEA "/home2/mosab/2017.01.04_V_MO/2017_12_12_V_STRAIN/annealing"	// restart file name
#define READ_FILEB "/home2/mosab/2017.01.04_V_MO/2017_12_12_V_STRAIN/annealing"	// restart file name

#define POT_FILE1 "/home/mosab/resources/V_BN_2017_mosab_test.eam.fs"		// potential file name   //// note this one changed 14/05/2017 
#define ATOM_SYMBOL1 "V"	// atomic symbol
#define DISCARD_NRUN1 0		// the number of run to change the timing of recoil event


void func_write_execution1_command(char proj1[],int discard_nrun);
void func_create_input_files(int i,char fold1[], char proj1[], int discard_nrun);
void func_make_data_file(int i,int SA, int SB, int SC,double LATC, char OUTF[]);
void func_read_execution1_command(int i,char fold2[],char proj1[],int discard_nrun);//,int SA);//,int SB,int SC,int LATC);


int main(void)
{
	int i,j,n,s,t,SA,SB,SC;
	int ns;
	int discard_nrun;
	char cline1[500];
	char fold1[300],fold2[300];
	char command1[300],command2[300],command3[300],OUTF[500];
	char proj1[300];
	char fname1[300];
	double LATC;
	FILE *fr,*fw;

	for(i=1;i<2;i++)
	{	LATC=3.0399;
		func_make_data_file(i,SA,SB,SC,LATC,OUTF);
	
	
   	for(t=1;t<(N_STEP+1);t++)
       	{	

         	discard_nrun=t*D_STEP;
              	printf("#STEP:\t%d\n",discard_nrun);
		sprintf(fold1,"STEP_%.2d",discard_nrun);
		sprintf(command1,"mkdir %s",fold1);
		system(command1);
		//sprintf(command2,"mkdir %s",fold2);
		//system(command2);
		sprintf(command3,"mkdir results");
		system(command3);

		
				sprintf(proj1,"step_%.2d",t);
				func_create_input_files(i,fold1, proj1, discard_nrun);
				func_write_execution1_command(proj1,discard_nrun);
				printf("Write excution code completed");
				printf("echo \"%s\"\n",proj1);
    			 	printf("date\n");
			
				func_read_execution1_command(i,fold2,proj1,discard_nrun);//,SA,SB,SC,LATC);
	}
		}


	return discard_nrun;
}


void func_create_input_files(int i,char fold1[], char proj1[], int discard_nrun)
{
	int p,j,k,SA,SB,SC;
	char cline1[500],fname1[300],READ_FILE1[500];
	FILE *fw,*fr;
	SA=i*A;
	SB=i*B;
	SC=i*C;

    	//sprintf(fname1,"%s/input.%s",fold1,proj1);
	sprintf(fname1,"input.%s",proj1);
	sprintf(READ_FILE1,"%s/data.fe_bcc%d%d%d",READ_FILEA,SA,SB,SC);

      	if( (fw=fopen(fname1,"w"))==NULL )      {  printf("error in open %s file\n",fname1);  exit(1);  }
     	fprintf(fw,"#Step 1 for relaxation");

   	if( (fr=fopen(INP_1_LAMMPS_REF,"r"))==NULL )      {  printf("error in open INP_1_LAMMPS_REF file\n");  exit(1);  }

      	for(p=0;p<NLINE_INP_1_LAMMPS_REF;p++)
      	{
          	fgets(cline1,500,fr);
              	if(p==1)        fprintf(fw,"variable read_file1 string    %s\n",READ_FILE1);
              	else if(p==2)   fprintf(fw,"variable pot_file1 string    %s\n",POT_FILE1);
              	else if(p==3)   fprintf(fw,"variable atom_symbol1 string %s\n",ATOM_SYMBOL1);
             	else if(p==4)  fprintf(fw,"variable temp equal         %d\n",TEMP);
                else    fprintf(fw,"%s",cline1);
      	} 
     	fclose(fr);
     	fclose(fw);
}


void func_write_execution1_command(char proj1[],int discard_nrun)
{
	int n;
	char command4[500],command5[300],command6[300],command7[300];

	printf("echo \"%s\"\n",proj1);
     	printf("date\n");

	//sprintf(command4,"cd STEP_01/");//,discard_nrun);
	//system(command4);

      	sprintf(command5,"%s <input.step_01 >output.voronoi_%s",LAMMPS_EXE,proj1);
	system(command5);
      //	printf(command6,"mv output.voronoi_%s result\n",proj1);
	//system(command6);

        sprintf(command7,"rm input.step_01");
	system(command7);

}

void func_read_execution1_command(int i,char fold2[],char proj1[],int discard_nrun)//,int SA,int SB,int SC,int LATC)

{
	int n,j,k,index,SA,SB,SC,m;
	double lattice_cons1[1001],lattice_cons2,LATC,lattice_cons[1001],total;
	double volume,side,LC;
	char READ_FILE2[500],fname2[300],fname1[300];
	char cline1[500],OUTF[300],command9[500];
	char fold1[300],command8[500];//fold2[300];
	SA=i*A;
	SB=i*B;
	SC=i*C;
	//lattice_cons[1001]=0.0;
	//volume=0.0;
	total=0;
	double po=1.0/3.0;
	FILE *fw1,*fw4,*fr1,*fr2;

    	sprintf(fname1,"%s/input.%s",fold1,proj1);
   	if( (fr1=fopen(READ_OUTPUT,"r"))==NULL )      {  printf("error in open READ_OUTPUT file\n");  exit(1);  }
	

      //	for(n=0;n<N_LINE_OUTPUT_FILE;n++)
	      //	{
	for(n=0;n<21;n++) {fgets(cline1,500,fr1);}
		
        	for(n=0;n<1001;n++)
		{fscanf(fr1,"%lf",&volume);
		//printf("volume %lf\n",volume);
		side=volume/(SA*SB*SC);

		lattice_cons[n]=pow(side,1.0/3.0);
 		total+=lattice_cons[n];
		printf("lattice  %.9f \n",lattice_cons[n]);
      		} 

	//fclose(fr1);
	LATC=total/1001;
		printf("SA=%d  SB=%d  SC=%d\n",SA,SB,SC);
	     	printf("lattice contstat=%lf\n",LATC);
	sprintf(READ_FILE2,"%s/data.fe_bcc%d%d%d",READ_FILEB,SA,SB,SC);

	func_make_data_file(i,SA, SB, SC,LATC, OUTF);  ///generate again data file using new LATC
	

	sprintf(fname2,"input.2");
      	if( (fw1=fopen(fname2,"w"))==NULL )      {  printf("error in open %s file\n",fname2);  exit(1);  }
     	fprintf(fw1,"#Step 2 for relaxation");

   	if( (fr2=fopen(INP_2_LAMMPS_REF,"r"))==NULL )      {  printf("error in open INP_2_LAMMPS_REF file\n");  exit(1);  }

      	for(n=0;n<NLINE_INP_2_LAMMPS_REF;n++)
      	{
          	fgets(cline1,500,fr2);
              	if(n==1)        fprintf(fw1,"variable read_file2 string    %s\n",READ_FILE2);
              	else if(n==2)   fprintf(fw1,"variable pot_file1 string    %s\n",POT_FILE1);
              	else if(n==3)   fprintf(fw1,"variable atom_symbol1 string %s\n",ATOM_SYMBOL1);
             	else if(n==4)  fprintf(fw1,"variable temp equal         %d\n",TEMP);
                else    fprintf(fw1,"%s",cline1);
      	} 
      	

     	fclose(fr2);
     	fclose(fw1);
	sprintf(command8,"%s <input.2 >output2.voronoi_%s",LAMMPS_EXE,proj1);
	system(command8);
	sprintf(command9,"cp restart.fe- restart.%d%d%d ",SA,SB,SC);
	system(command9);

    
}





void func_make_data_file(int i,int SA, int SB, int SC,double LATC, char OUTF[])
{
        //int i;
        int x,y,z,r;//SA,SB,SC;
	//char OUTF[500];
	//double LATC=3.16;	
        const double rpos[NR][3]={
                {0.0, 0.0, 0.0},
                {0.5, 0.5, 0.5}
        };
	SA=A*i;
	SB=B*i;
	SC=C*i;
	printf("SA=%d   SB=%d   SC=%d",SA,SB,SC); 
	sprintf(OUTF,"data.fe_bcc%d%d%d",SA,SB,SC);
        FILE *fw3;


        if( (fw3=fopen(OUTF, "w"))==NULL )
        {
                printf("error: cannot open the file\n");
                exit(1);
        }


        fprintf(fw3, "# Fe Structure bcc %dx%dx%d #\n\n", SA,SB,SC);
        fprintf(fw3, "%d   atoms\n", NR*SA*SB*SC);
        fprintf(fw3, "1  atom types\n\n");
        fprintf(fw3, "%.10f  %.10f   xlo xhi\n", 0.00, LATC*SA);
        fprintf(fw3, "%.10f  %.10f   ylo yhi\n", 0.00, LATC*SB);
        fprintf(fw3, "%.10f  %.10f   zlo zhi\n\n", 0.00, LATC*SC);
        fprintf(fw3, "Masses\n\n");
        fprintf(fw3, "1  %.10lf\n\n",MASS);
        fprintf(fw3, "Atoms\n\n");

        i=0;
        for(x=0;x<SA;x++)
        {
                for(y=0;y<SB;y++)
                {
                        for(z=0;z<SC;z++)
                        {
                                for(r=0;r<NR;r++)
                                {
                                        i++;
                                        fprintf(fw3, "%d  %d  %.12lf  %.12lf  %.12lf\n",
                                                 i, 1, (x+rpos[r][0])*LATC, (y+rpos[r][1])*LATC, (z+rpos[r][2])*LATC);
                                }
                        }
                }
        }

        if(i!=(SA*SB*SC*NR))    {  printf("error in the number of atoms\n");  exit(1);  }

        fclose(fw3);

}








