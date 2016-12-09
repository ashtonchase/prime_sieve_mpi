//MIT License

//Copyright (c) 2016 Ashton Johnson

//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:

//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.

//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

/******************************************************************/
/**
* Prime number generation program              -- MPI Daisy-Chain version 
* @file prime_mpi.cpp
* @author Ashton Johnson, Paul Henny
* @date 15 October
* @brief MPI Version of Prime number generation program 
* This is an implementation utilizes collective communication to 
* distribute prime number of which multiples should be excluded  
* from the prime number list.                                      
* 
*/
// prime.cpp
// compilation:
//   gnu compiler
//      g++ prime.cpp -o prime -O3 -lm
// Note: to compile a parallel MPI program version which is named
//   prime_mpi.cpp
//   then execute the following command
//      gnu compiler
//         mpic++ prime_mpi.cpp -o prime_mpi -lm  -O3
/*
  Parallel data passing model.  The prime to check the data sample
  is passed to each process after the previous finishes.

  To execute:
  prime_mpi max_numb
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
#include <sys/time.h>
#include <mpi.h>


/*! PRIME_EXIT is value passed to indicated there are not more values to 
  / distribute from the master to slave processes. 
*/
#define PRIME_EXIT (int)-1

/*! MPI Tag used for communication
 */
#define MARK_PRIME_TAG (int)7

/*! Timer related definitions 
 */
#define TIMER_CLEAR     (tv1.tv_sec = tv1.tv_usec = tv2.tv_sec = tv2.tv_usec = 0)
#define TIMER_START     gettimeofday(&tv1, (struct timezone*)0)
#define TIMER_ELAPSED   ((tv2.tv_usec-tv1.tv_usec)+((tv2.tv_sec-tv1.tv_sec)*1000000))
#define TIMER_STOP      gettimeofday(&tv2, (struct timezone*)0)
struct timeval tv1,tv2;

/// highest number passed from the user to check all numbers for prime eligibility 
int highestNumber;
/// square root of highestNumber
int rootHighestNumber;
/// pointer for value of the local array size
int *localArraySize;
/// pointers for arrays of all number and the local processes array
int *numberArray, *localNumberArray;
/// pointers for arrays indicating whether or not the number is prime
bool *isPrimeArray, *lclIsPrimeArray;
/// MPI Specifics for the number or processes and the rank
int numProc, myRank;
MPI_Comm   *mpiPrimeComm;
MPI_Group  *world_group;

/**
   Routine to retrieve the highest number to search for all lower valued possibilities of prime numbers
*/
void GetMaxNumber(int argc,char *argv[],int *highestNumber) {
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

/**
  Routine that outputs the matrices to the screen 
*/
void PrintMatrix(float *array,int dim_m,int dim_n)
{
  int i,j;
  for(i=0;i<dim_m;i++) {
    for (j=0;j<dim_n;j++) {
      cout << array[i*dim_n+j] << " ";
    }
    cout << endl;
  }
}

/** \brief Sieves through local numbers to mark off primes.
 * \param myRank MPI rank of the local process within. 
 * \param numProc MPI total number of proccesses.
 * \param localArraySize amount of numbers covered by the local array.
 * \param isPrimeArray boolean array for the numbers covers in the 
 * local process
 * 
 * Detailed Description
 */
void ComputePrimes(int myRank, int numProc,  int *localArraySize, bool  isPrimeArray[])
{
#ifdef DEBUG
  cout<<"Rank:"<<myRank<<"\tComputing Primes."<<endl;
#endif     
  /// Holds the recieved value of mutliples to mark off
  int *primeToMark= new (nothrow) int;
  /// Dynamic lowest index of the local array. Increases as lowest value
  /// is marked off
  int lowestIndex=myRank*(*localArraySize);
  /// Initial, unchanging lowest index of the local array
  const int baseIndex=lowestIndex;
  /// Dynamic highest index of the local array. Increases as lowest values
  int highestIndex=myRank*(*localArraySize)+(*localArraySize);
  /// Initial, unchangin highest index of the local array
  const int baseUpperIndex=highestIndex;
  /// lowestMultiple is used to set lowest bounds of for loop that will be
  /// used to mark off values
  int lowestMultiple=lowestIndex;
  /// highestMultiple is used to set lowest bounds of for loop that will be
  /// used to mark off values
  int highestMultiple=highestIndex;

  if (myRank==0){
    /// Set highest index. Initially it will be the root highest number,
    /// but that might be lowered based upon what gets marked off.
    highestIndex=min(highestIndex,rootHighestNumber);
    /// Start sending numbers from Rank 0.
    for ( int i=2; i<highestIndex; i++){
#ifdef DEBUG
      cout<<endl<<"Checking: "<<i;
#endif
      /// Check to see if the current number is a prime. 
      if (isPrimeArray[i]==true){
#ifdef DEBUG
	cout<<" Sending..."<<endl;
#endif
	/// Send number to all other ranks.
	*primeToMark=i;
	for( int dst=myRank+1; dst<numProc;dst++)
	  MPI_Send(primeToMark, 1, MPI_INT, dst, MARK_PRIME_TAG, MPI_COMM_WORLD);
    
	/// Update local primes.
	for (int j=i*2 ; j<baseUpperIndex ; j+=i)
	  isPrimeArray[j]=false;   //mark all multiples as non-primes    
      
      }
    }
#ifdef DEBUG
    cout<<"Sending EXIT..."<<endl;
#endif
    /// Send exit by sending the PRIME_EXIT constat to other proccesses.
    *primeToMark=PRIME_EXIT;
    for( int dst=myRank+1; dst<numProc;dst++) 
      MPI_Send(primeToMark, 1, MPI_INT, dst, MARK_PRIME_TAG, MPI_COMM_WORLD);
#ifdef DEBUG
    for ( int i=2; i<baseUpperIndex; i++){    
      if (isPrimeArray[i-baseIndex]==true ) 
	cout<<"R:"<<myRank<<"\t"<<i<<" is Prime."<<endl;
    }
#endif
    return;
    
    
    
    /// Slave Thread Work
  }else{   
    /// Run indefinitely until receiving PRIME_EXIT 
    while(1)
      {
    
#ifdef DEBUG
	for ( int i=lowestIndex; i<highestIndex; i++)
	  cout<<i<<"\t"<<(int)isPrimeArray[i-baseIndex]<<endl;
	cout<<"Rank:"<<myRank<<"\tloMultiple:"<<lowestMultiple<<"\thiMultiple:"<<highestMultiple;
	cout<<"\tloIndex:"<<lowestIndex;
	cout<<"\thiIndex:"<<highestIndex<<endl;
	cout<<"R:"<<myRank<<"  Waiting to Rx..."<<endl;
#endif
	/// Blocking receive from master process.
	MPI_Recv(primeToMark,1, MPI_INT, 0, MARK_PRIME_TAG, MPI_COMM_WORLD,NULL);
      
#ifdef DEBUG
	cout<<"R:"<<myRank<<" Received: "<<*primeToMark<<endl;
#endif
        /// Check to see if valid number or exit flag
        if(*primeToMark==PRIME_EXIT){
#ifdef DEBUG
	  cout<<"R:"<<myRank<<"EXIT Received"<<endl;        
	  for ( int i=baseIndex; i<baseUpperIndex; i++){    
	    if (isPrimeArray[i-baseIndex]==true )
	      cout<<"R:"<<myRank<<"\t"<<i<<" is Prime."<<endl;
	  }
#endif     
	  /// Exit while loop
	  return;
	}     

	/// Determine lowest and highest multiple of the primeToMark sent 
	/// that is in our array.
	if ( lowestIndex % (*primeToMark) !=0 ){
	  lowestMultiple=(lowestIndex/(*primeToMark)+1)*(*primeToMark);
	}else{
	  lowestMultiple=lowestIndex;
	} 

	/// Trim up the lowestIndex as required until one that hasn't
	/// already been marked is found
	for ( int z=lowestIndex; z<highestIndex; z++){    
	  if (isPrimeArray[z-baseIndex]==false ){
	    lowestIndex++;
#ifdef DEBUG
	    cout<<"incrementing lowestIndex"<< lowestIndex;        
	    cout<<" because isPrimeArray[i-lowestIndex]==false=="<<(int)isPrimeArray[z-baseIndex]<<endl;
#endif
	  } else { z=highestIndex;}
	}
       
     

	/// Mark all multiples as non-primes    
	for ( int i=lowestMultiple; i<highestIndex; i+=(*primeToMark)){
	  /// If we just marked off the lowest index, reset the lowest index 
	  isPrimeArray[i-baseIndex]=false;  
#ifdef DEBUG
	  cout<<"Setting "<<i<<" to false"<<endl;
#endif
         
	}
	/// Trim down the highest Index as required
	for ( int i=(highestIndex-1); i>lowestIndex; i--){    
	  if (isPrimeArray[i-baseIndex]==false ){ highestIndex--;} else { break;}
	}

      }//while(1)
  }
}//compute


/**
   \param argc input argument characater count.
   \param argv input argument character array.
*/
int main( int argc, char *argv[])
{
  
  /// MPI communication type pointer
  MPI_Comm *mpiPrimeComm_t;
  /// Initialize MPI environment
  MPI_Init(&argc,&argv);
  /// Get total number of MPI proceesses
  MPI_Comm_size(MPI_COMM_WORLD,&numProc); 
#ifdef DEBUG
  cout<<"Number of Processes: "<<numProc<<endl;
#endif
  /// Get unique task id number    
  MPI_Comm_rank(MPI_COMM_WORLD,&myRank); 
#ifdef DEBUG
  cout<<"Processes Rank: "<<myRank<<endl;
#endif    

#ifdef DEBUG
  cout<<"MPI Comm Created"<<endl;
#endif    
  /// Initialize local array size
  localArraySize= new (nothrow) int;
  if (localArraySize==0) exit(1);


  
  /// Get matrix sizes
  GetMaxNumber(argc,argv,&highestNumber);

  /// Determine square root
  rootHighestNumber=sqrt(highestNumber);
#ifdef DEBUG
  cout<<"Rank:"<<myRank<<"\tRoot of Highest Number: "<<rootHighestNumber<<endl;
#endif  
  /// Determine the local group size
  if (myRank<(highestNumber%numProc)) {
    *localArraySize=(highestNumber/numProc)+1;
  }else{ 
    *localArraySize=(highestNumber/numProc);
  }
#ifdef DEBUG
  cout<<"Rank:"<<myRank<<"\tLocal Array Size:"<<*localArraySize<<endl;
#endif    
  
  /// Allocate the local prime array
  lclIsPrimeArray = new (nothrow) bool[*localArraySize];
  /// Test for correct allocation
  if(lclIsPrimeArray==0) {
    cout <<"Rank:"<<myRank<<"\tERROR:  Insufficient Memory" << endl;
    exit(1);
  }
#ifdef DEBUG
  cout<<"Rank:"<<myRank<<"\tlclIsPrimeArray @ : 0x"<<hex<<lclIsPrimeArray<<dec<<endl;
#endif        
  /// Initialize isPrimeArray to true, all non-primes
  /// will be marked false.
  for ( int i=0;i<*localArraySize;i++) lclIsPrimeArray[i]=true;
#ifdef DEBUG
  cout<<"Rank:"<<myRank<<"\tLocal Prime Array Initialized."<<endl;
#endif         
 



  /***--------Compute Primes, Locally------------
   *In this routine, the lowest rank will start sending out numbers in groups of 1
   *that will be prime numbers that subsequent ranks can mark off.
   *
   *If rank contains numbers that are less than the square root of the highest number
   *then we will need to do work, sequentially, first letting rank zero complete and
   *signal the next task to begin sending out prime numbers to subsequent processes
   *above it. 
   */
  
  TIMER_CLEAR;    
  TIMER_START;
  /// Call function on all processes to seive through the primes.
  ComputePrimes(myRank,numProc,localArraySize, lclIsPrimeArray);
  MPI_Barrier(MPI_COMM_WORLD);
  
  TIMER_STOP;  

#ifdef PRINT_PRIMES
  /// Print primes
  for ( int i=1;i<highestNumber;i++) 
    if (isPrimeArray[i]==true)cout<<i<<endl;
#endif
 
  if (myRank==0)
    cout << "time=" << setprecision(8) <<  TIMER_ELAPSED/1000000.0  << " seconds" << endl;
  /// Terminate MPI communications
  MPI_Finalize();
}


