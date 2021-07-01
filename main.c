#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <time.h>
#include <memory.h>
#define nmb_of_cities 8
#define DATAFILE "cities.txt"


typedef struct 
{
    float lon;
    float lat;
} city;

int POPULATION_SIZE;  //=1000;
int GENERATIONS;      // =25;
float MUTATIONFACTOR; // =0.04;
int EXCHANGEPERIOD;   // =10;

float globalfit;

MPI_Datatype MPI_City;
typedef int path[nmb_of_cities];
city cities[nmb_of_cities];
void server(void);
void slave(void);
int rank, size;

// local_fittest commmand-line arguments order
// pop size, generations, mutation factor, exchangeperiod



void shuffle_cities_in_path(path p)
{
    int i;
    int a, b;
    int temp;
    for (i = 0; i < nmb_of_cities; i++)
    {
        a = rand() % nmb_of_cities;
        b = rand() % nmb_of_cities;
        temp = p[a];
        p[a] = p[b];
        p[b] = temp;
    }
}

float square(float x)
{
    return x * x;
}
float calculate_distance(city c1, city c2)
{
    return sqrt(square(c1.lon - c2.lon) +
                square(c1.lat - c2.lat));
}
float calculate_fitness(path p)
{
    int i;
    float sum = 0;
    for (i = 0; i < nmb_of_cities; i++)
        sum += calculate_distance(cities[p[i]], cities[p[(i + 1) % nmb_of_cities]]);
    return sum;
}
int comparepaths(const void *p1, const void *p2)
{
    float f1, f2;
    f1 = calculate_fitness(*(path *)p1);
    f2 = calculate_fitness(*(path *)p2);
    if (f1 < f2)
        return 1;
    else
        return (f1 == f2) ? 0 : -1;
}
void crossover(path p1, path p2)
{
    int cp, i, j, k;
    path newp1, newp2;
    cp = rand() % nmb_of_cities;
    for (i = 0; i < cp; i++)
        newp1[i] = p1[i];
    //start p2 from the beginning
    // if an element does not exist add it
    j = 0;
    while (i < nmb_of_cities)
    {
        for (k = 0; k < i; k++)
            if (p2[j] == newp1[k])
                break;
        if (k == i)
            newp1[i++] = p2[j];
        j++;
    }
    for (i = 0; i < cp; i++)
        newp2[i] = p2[i];
    //start p1 from the beginning
    // if an element does not exist add it
    j = 0;
    while (i < nmb_of_cities)
    {
        for (k = 0; k < i; k++)
            if (p1[j] == newp2[k])
                break;
        if (k == i)
            newp2[i++] = p1[j];
        j++;
    }
    memcpy(p1, newp1, sizeof(newp1));
    memcpy(p2, newp2, sizeof(newp2));
}

void display_path(path p)
{
    int i, f;
    for (i = 0, f = 1; i < nmb_of_cities; i++)
        f *= (p[i] == 0) ? 1 : p[i];
    for (i = 0; i < nmb_of_cities; i++)
        printf("%d   ", p[i]);
    printf("\nfitness= %f\n", calculate_fitness(p));
}

void mutate(path p)
{
    int c1, c2;
    int c;
    c1 = rand() % nmb_of_cities;
    c2 = rand() % nmb_of_cities;
    c = p[c1];
    p[c1] = p[c2];
    p[c2] = c;
}
void select_crossover_choices(int *p, float *fits)
{
    int i, j;
    float placeholder, sum = 0;
    placeholder = fits[POPULATION_SIZE - 1];
    //implements the selection operator
    //pre: fits array is sorted largest first. less is better
    for (i = POPULATION_SIZE - 2; i >= 0; i--)
    {
        fits[i] += fits[i + 1];
    }
    for (i = 0; i < 2 * (POPULATION_SIZE - 2); i++)
    {
        //get candidates for crossover
        //making it a bigger probability to crossover those with longer distances (less needed) and less probability to crossover less distances
        float r = rand() / RAND_MAX;
        for (j = 0; j < POPULATION_SIZE; j++)
            if (fits[j] < r)
                break;
        p[i] = j - 1;
    }
}



int main(int argc, char **argv)
{
    int i;
   
    city city_one;
    int len_struct_var[2];
    MPI_Aint loc_struct_var[2];
    MPI_Datatype type_struct_var[2];
    MPI_Aint baseaddr;

     MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     MPI_Comm_size(MPI_COMM_WORLD, &size);

    //start timer
    double start_time = MPI_Wtime();

    //get variables
    POPULATION_SIZE = atoi(argv[1]);
    GENERATIONS = atoi(argv[2]);
    MUTATIONFACTOR = atof(argv[3]);
    EXCHANGEPERIOD = atoi(argv[4]);
    // Make the city MPI type
    MPI_Get_address(&city_one, &baseaddr);
    MPI_Get_address(&city_one.lon, &loc_struct_var[0]);
    MPI_Get_address(&city_one.lat, &loc_struct_var[1]);

    len_struct_var[0] = 1;
    loc_struct_var[0] -= baseaddr;
    type_struct_var[0] = MPI_FLOAT;
    len_struct_var[1] = 1;
    loc_struct_var[1] -= baseaddr;
    type_struct_var[1] = MPI_FLOAT;

    //creating the MPI_STRUCT with previously defined parameters
    MPI_Type_create_struct(2, len_struct_var, loc_struct_var, type_struct_var, &MPI_City);
    MPI_Type_commit(&MPI_City);
    if (rank == 0)
    {
        int i;
        FILE *fcities;
    fcities = fopen(DATAFILE, "r");
      
        for (i = 0; i < nmb_of_cities; i++)
        {
            fscanf(fcities, "%f %f", &cities[i].lon,
                   &cities[i].lat);
        }
        fclose(fcities);
        // with nothing else to do, the server can work as a slave
        slave();
    }
    else

        slave();

    // Finish up
    if (rank == 0)
    {
        double end_time = MPI_Wtime();
        printf("Wallclock time elapsed: %.8lf seconds\n", end_time - start_time);
    }

     MPI_Finalize();
    return 0;
}



void slave(void)
{
    int i, j;
    path *paths, *placeholder;
    float *pathfits;
    

    int crossover_choices[2 * (POPULATION_SIZE - 2)];
    MPI_Status status;
    paths = (path *)malloc(POPULATION_SIZE * sizeof(path));
    pathfits = (float *)malloc(POPULATION_SIZE * sizeof(float));
    MPI_Bcast(cities, nmb_of_cities, MPI_City, 0, MPI_COMM_WORLD);
    // Initialize random number generator
    
    srand(time(0));
    // -------------------------------------- generate initial paths ----------------------------------------------------------------------------

    for (i = 0; i < nmb_of_cities; i++)
        paths[0][i] = i;
    for (i = 1; i < POPULATION_SIZE; i++)
    {
        memcpy(paths[i], paths[i - 1], sizeof(paths[i]));
        shuffle_cities_in_path(paths[i]);
    }
    // -------------------------------------- looping number of generations times --------------------------------------
    for (j = 0;; j++)
    {
        // sort paths (the fittest paths will be local_fittest the end)
        qsort(*paths, POPULATION_SIZE, sizeof(paths[0]), comparepaths);
        if (j == GENERATIONS)
            break;

        // so that no process after many generations is left wayy behind and is kept up with who is the fittest path and places it local_fittest the place of the least fittest path
        if (j % EXCHANGEPERIOD == 0)
            MPI_Sendrecv(paths[POPULATION_SIZE - 1], nmb_of_cities, MPI_INT,
                         (rank + 1) % size, 0, paths[0], nmb_of_cities, MPI_INT,
                         (rank - 1) % size, 0, MPI_COMM_WORLD, &status);

        //calculate the fitness of every path and place it local_fittest its index local_fittest pathfits
        for (i = 0; i < POPULATION_SIZE; i++)
            pathfits[i] = calculate_fitness(paths[i]);
        select_crossover_choices(crossover_choices, pathfits);
      

        for (i = 0; i < POPULATION_SIZE - 2; i += 2)
            crossover(paths[crossover_choices[i]],
                      paths[crossover_choices[i + 1]]);
        for (i = 0; i < POPULATION_SIZE - 2; i++)
            if (rand() / RAND_MAX < MUTATIONFACTOR)
                mutate(paths[i]);
    }
    
    // printf("\nprocess %d found Best calculate_fitness is:%f\n",rank,calculate_fitness(paths[POPULATION_SIZE - 1]));
    float localfit = calculate_fitness(paths[POPULATION_SIZE - 1]);
    // MPI_Allreduce(&local_fittest, &global_fittest, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    MPI_Reduce(&localfit, &globalfit, 1, MPI_FLOAT, MPI_MIN,0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("process %d found Best calculate_fitness is: %f\n",rank,globalfit);
        //display_path(paths[POPULATION_SIZE - 1]);
    }
    free(pathfits);
    free(paths);
}
