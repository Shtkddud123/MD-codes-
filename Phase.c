#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

int main ()
{


 /*box dimensions*/

  double box1;
  double box2; 
  double boxlength[boxdim];


  /*Each line from the dump file*/
  char line[100];
  int i = 0;
  int numberofatoms = 17828; 
  int atomtype;


  /*Misc Parameters*/ 

  int index; 
  int i = 0;
  int atomtype;


  double xco[numberofatoms],yco[numberofatoms],zco[numberofatoms]; 
  int a[numberofatoms],b[numberofatoms];
  double x,y,z; /*coordinates for the atoms in the box*/

  //we want to allocate enough memory so that we can count in all the possible lines inside a dump file. 
  /*open the file for reading*/
  FILE *NewFile;
  NewFile = fopen("/home/chem/msrgbj/cProgrammingCodes/dump.c12e2_Nchain500_Nwater14328_npt", "r");


  //The total number of lines we have is 3585238
  //Need to change this into a malloc as we will be looking at lots of dumpfiles with different lengths

  //going through all lines i 
  for (i=0;i<numberoflines;i++) {


    // int l looks for all the strings 

    int l = 0; 

    if(NewFile == NULL) {

      printf("Error opening file\n");
      exit(1);
    }

    fgets(line,sizeof(line),NewFile);

    if (l < 5) {

      //acknowledging the parts of the dump file where we know there is no useful information 

    }

    else if ((l > 4 && l < 8)) {
      
	/*We are scanning the bit with just the box parameters*/
      
	sscanf(line, "%lf %lf", &box1, &box2);
	//	printf("%lf %lf \n",box1,box2 );
	boxlength[l-5] = box2-box1; 
	l++;
	
    }

    else if (l == 8) {
      
	//	printf(" **** l= %d \n",l);
	/*We are doing nothing*/
      
	l++;
      }
    

    else { 

      sscanf(line,"%d %d %lf %lf %lf",&index,&atomtype,&x,&y,&z);
	a[index-1] = index;
	b[index-1] = atomtype;
	xco[index-1] = x*boxlength[0]; 
	yco[index-1] = y*boxlength[1];
	zco[index-1] = z*boxlength[2];  
	n++;
	l++;

    }




    printf("%s ",line);

  }



}
