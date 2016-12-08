/******************************************************************/
/* Prime number generation program              -- MPI Daisy-Chain version */
/* 6 November 2016 */
/*Copyright 2016 Ashton Johnson, Paul Henny */
/******************************************************************/
// prime.cpp
// compilation:
//   gnu compiler
//      g++ prime.cpp -o prime -O3 -lm
// Note: to compile a parallel MPI program version which is named
//   prime_chain.cpp
//   then execute the following command
//      gnu compiler
//         mpic++ prime_chain.cpp -o prime_chain -lm  -O3
/*
  Daisy-chain data passing model.  The prime to check the data sample
  is passed to each process after the previous finishes.

  To execute:
  prime_chain max_numb
*/

//#define TESTING
using namespace std;
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h> // for MPI parrallelism
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
  Routine to retrieve the highest number to search for all lower valued possibilites of prime numbers
*/
void get_max_number(int argc,char *argv[],int *highestNumber) {
  if(argc!=2) {//if it is not two arguments (including the program name)
    cout<<"usage:  prime <highestNumber>"
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
  MAIN ROUTINE: summation of a number list
*/

int main( int argc, char *argv[])
{

  int highestNumber;
  int rootHighestNumber;
  int *numberArray;
  bool *isPrimeArray;
  MPI_Status status;
  int numtasks,rank, num_to_send;
  int curr_prime = 2;
  int last_nonprime;

  MPI_Init(&argc,&argv); // initialize MPI environment
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks); // get total number of MPI processes
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); // get unique rank of the process

  int rec_prime, rec_lastnon = 2;
  int type;
  bool *prime_buf;

  rec_prime = 2;

  /*
     get matrix sizes
  */
  get_max_number(argc,argv,&highestNumber);

  // The amount to send to each process
  num_to_send = ceil((double)(highestNumber+1)/(double)numtasks);
  //cout << highestNumber << endl;
  //cout << numtasks << endl;

  //determine square root
  rootHighestNumber=ceil(sqrt(highestNumber+1));
  //cout << "the highest nubmer "<< rootHighestNumber << endl;

  // dynamically allocate
  prime_buf    = new (nothrow) bool[num_to_send];
  // test for correct allocation
  if(prime_buf==0) {
    cout <<"ERROR:  Insufficient Memory" << endl;
    MPI_Finalize(); // Exit MPI
    exit(1);
  }

  if(rank==0) {
    // dynamically allocate
    //isPrimeArray = new (nothrow) bool[(highestNumber+1)+(highestNumber+1)%numtasks];
    isPrimeArray = new (nothrow) bool[num_to_send*numtasks];
    // test for correct allocation
    if(isPrimeArray==0) {
      cout <<"ERROR:  Insufficient Memory" << endl;
      MPI_Finalize(); // Exit MPI
      exit(1);
    }
    //initialize isPrimeArray to true, all non-primes
    //will be marked false.
    for ( int i=0;i<highestNumber+1;i++) isPrimeArray[i]=true;
  }

  /*
    Start recording the execution time
  */
  TIMER_CLEAR;
  TIMER_START;

  //This outer loop will go through all the numbers up
  //to the sqaure root of the highest number chosen.
  //At each number, if the number has not been marks as
  //non-prime in the previous iterations, it is prime.


  /// Scatter isPrimeArray, put 0 padding on end
  MPI_Scatter(isPrimeArray,num_to_send,MPI::BOOL,prime_buf,num_to_send,MPI::BOOL,0,MPI_COMM_WORLD);

  /// Broadcast the size of numbers sent to each
  MPI_Bcast(&num_to_send,1,MPI_INT,0,MPI_COMM_WORLD);

  if (rank==0) {
    cout << "Largest Number:  " << highestNumber << endl; // Debug
    type = 123;
    for(int i=2;i<rootHighestNumber;i++) {
      curr_prime = i;
      last_nonprime = i;
      //check to see if the current number is a prime.
      if (prime_buf[i]==true){
        //mark all multiples as non-primes
        for (int j=i*2 ; j<num_to_send ; j=j+i){
  	       prime_buf[j]=false;
           last_nonprime = j;
        }
        if (numtasks>1) {
          // Send what index that was just done to next process, and teh last non prime number
          MPI_Send(&curr_prime,1,MPI_INT,rank+1,type,MPI_COMM_WORLD);
          MPI_Send(&last_nonprime,1,MPI_INT,rank+1,type,MPI_COMM_WORLD);
        }
      }
    }
    curr_prime = rootHighestNumber+1;
    last_nonprime = rootHighestNumber+1;
    
    if (numtasks>1) {
      // Send a large number to indicate that we have done the entire array
      // This indicates that we are done TESTING
      // TODO: handle if the last prime is not in the first process's numbers
      MPI_Send(&curr_prime,1,MPI_INT,rank+1,type,MPI_COMM_WORLD);
      MPI_Send(&last_nonprime,1,MPI_INT,rank+1,type,MPI_COMM_WORLD);
    }
  }
  else {
    while(rec_prime<rootHighestNumber) {
      type = 123 + rank - 1;
      MPI_Recv(&rec_prime,1,MPI_INT,rank-1,type, MPI_COMM_WORLD,&status);
      MPI_Recv(&rec_lastnon,1,MPI_INT,rank-1,type, MPI_COMM_WORLD,&status);

      // Only do array math if it is a valid number
      if (rec_prime<rootHighestNumber) {
        // Debug strings
        //cout << "rank " << rank << " Rec Prime " << rec_prime << endl;
        //cout << "rank " << rank << " Last Not " << rec_lastnon << endl;

        //mark all multiples as non-primes
        for(int j=rec_lastnon+rec_prime-(num_to_send*rank); j<num_to_send; j=j+rec_prime) {
           prime_buf[j]=false;
           rec_lastnon = j+num_to_send*rank;
        }
      }
      // Send the prime and last non-prime to next process, if it isn't the last
      if(rank<numtasks-1) {
        type = 123 + rank;
        // Send what index that was just done to next process
        MPI_Send(&rec_prime,1,MPI_INT,rank+1,type,MPI_COMM_WORLD);
        MPI_Send(&rec_lastnon,1,MPI_INT,rank+1,type,MPI_COMM_WORLD);
      }
    }
    // Print individual processes arrays
    //cout << "rank = " << rank << endl;
    //for ( int i=0;i<num_to_send;i++)
    //  if (prime_buf[i]==true)cout<<rank<<" ===== "<<i<<endl;
  }

   /// MPI Gather of primeArray
   MPI_Gather(prime_buf,num_to_send,MPI::BOOL,isPrimeArray,num_to_send,MPI::BOOL,0,MPI_COMM_WORLD);

  /*
    stop recording the execution time
  */
  TIMER_STOP;

  if(rank==0) {
    //print primes
    //for ( int i=1;i<highestNumber+1;i++)
    //  if (isPrimeArray[i]==true)cout<<i<<endl;

    cout << "time=" << setprecision(8) <<  TIMER_ELAPSED/1000000.0  << " seconds" << endl;
  }

  MPI_Finalize(); // Exit MPI
}
