/******************************************************************/
/* Prime number generation program              -- serial version */
/* 15 October 2016 
   /*Copyright 2016 Ashton Johnson, Paul Henny */
/******************************************************************/
// mm_mult_serial.cpp
// compilation:
//   gnu compiler
//      g++ prime.cpp -o prime -O3 -lm

//#define TESTING
using namespace std;
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>



#define MX_SZ 320
#define SEED 2397           /* random number seed */
#define MAX_VALUE  100.0    /* maximum size of array elements A, and B */

/* copied from mpbench */
#define TIMER_CLEAR     (tv1.tv_sec = tv1.tv_usec = tv2.tv_sec = tv2.tv_usec = 0)
#define TIMER_START     gettimeofday(&tv1, (struct timezone*)0)
//#define TIMER_ELAPSED   (tv2.tv_sec*1e6 + tv2.tv_usec) - (tv1.tv_sec*1e6 + tv1.tv_usec)
#define TIMER_ELAPSED   ((tv2.tv_usec-tv1.tv_usec)+((tv2.tv_sec-tv1.tv_sec)*1000000))
#define TIMER_STOP      gettimeofday(&tv2, (struct timezone*)0)
struct timeval tv1,tv2;

/*
  This declaration facilitates the creation of a two dimensional 
  dynamically allocated arrays (i.e. the lxm A array, the mxn B
  array, and the lxn C array).  It allows pointer arithmetic to 
  be applied to a single data stream that can be dynamically allocated.
  To address the element at row x, and column y you would use the
  following notation:  A(x,y),B(x,y), or C(x,y), respectively.
  Note that this differs from the normal C notation if A were a
  two dimensional array of A[x][y] but is still very descriptive
  of the data structure.

*/



/*
  Routine to retrieve the highest number to search for all lower valued possibilites of prime numbers
*/
void get_max_number(int argc,char *argv[],int *highestNumber) {
  if(argc!=2) {//if it is not two arguments (including the program name)
    cout<<"usage:  prime <highestNumber"
	<< endl;
    exit(1);
  }
  else {
    
    *highestNumber = atoi(argv[1]);
       
  }
  
  if (*highestNumber<=2 ) {
    cout<<"Error: highest number must be greater than 2"
	<< endl;
    exit(1);
  }
}

/*
  Routine that fills the number matrix with Random Data with values
  between 0 and MAX_VALUE
  This simulates in some way what might happen if there was a 
  single sequential data acquisition source such as a single file
*/
void fill_matrix(float *array,int dim_m,int dim_n)
{
  int i,j;
  for(i=0;i<dim_m;i++) {
    for (j=0;j<dim_n;j++) {
      array[i*dim_n+j]=drand48()*MAX_VALUE;
    }
  }
}

/*
  Routine that outputs the matrices to the screen 
*/
void print_matrix(float *array,int dim_m,int dim_n)
{
  int i,j;
  for(i=0;i<dim_m;i++) {
    for (j=0;j<dim_n;j++) {
      cout << array[i*dim_n+j] << " ";
    }
    cout << endl;
  }
}

/*
  MAIN ROUTINE: summation of a number list
*/

int main( int argc, char *argv[])
{

  int highestNumber;
  int rootHighestNumber;
  int *numberArray;
  bool *isPrimeArray;

  /* 
     get matrix sizes
  */
  get_max_number(argc,argv,&highestNumber);
  //determine square root
  rootHighestNumber=sqrt(highestNumber);



  // dynamically allocate 
  isPrimeArray = new (nothrow) bool[highestNumber];
  // test for correct allocation
  if(isPrimeArray==0) {
    cout <<"ERROR:  Insufficient Memory" << endl;
    exit(1);
  }
  //initialize isPrimeArray to true, all non-primes
  //will be marked false.
  for ( int i=0;i<highestNumber;i++) isPrimeArray[i]=true;

  /*
    Start recording the execution time
  */
  TIMER_CLEAR;
   TIMER_START;


  
  //This outer loop will go through all the numbers up
  //to the sqaure root of the highest number chosen.
  //At each number, if the number has not been marks as
  //non-prime in the previous iterations, it is prime.

  for ( int i=2; i<rootHighestNumber; i++){
    //check to see if the current number is a prime. 
    if (isPrimeArray[i]==true){
      //mark all multiples as non-primes
      for (int j=i*2 ; j<=highestNumber ; j=j+i){
	isPrimeArray[j]=false;
      }
    }

  }
    TIMER_STOP;  

  //print primes
    //  for ( int i=1;i<highestNumber;i++) 
    // if (isPrimeArray[i]==true)cout<<i<<endl;

 
  /*
    stop recording the execution time
  */ 

 
  cout << "time=" << setprecision(8) <<  TIMER_ELAPSED/1000000.0  << " seconds" << endl;
}


