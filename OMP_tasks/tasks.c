#include <stdio.h>
int main(int argc, char *argv[]) {
    #pragma omp parallel 
    { 
        #pragma omp single 
        {
            printf("A "); 
            #pragma omp task 
            {printf("race ");} 
             #pragma omp task 
            {printf("car ");} 
//            #pragma omp taskwait 
//            printf(“is fun to watch “);
        }
    } // end of parallel region
    printf("\n");
    return(0);
}
