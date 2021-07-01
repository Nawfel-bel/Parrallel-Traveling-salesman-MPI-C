#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void random_cities(int n, char outputFile[])
{
    int i;
    FILE *fid;
    fid = fopen(outputFile, "w");

    // the seed is the current time 
    srand(time(0));

     for( i = 0 ; i < n ; i++ ) {
      fprintf(fid, "%d", rand() % 300);
      fprintf(fid, " %d\n", rand() % 300);
   }
}


int main(int argc, char* argv[])
{   
    random_cities(30,"cities.txt");
    return 0;
}