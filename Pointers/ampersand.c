#include<stdio.h>

int* TAB() {
   int temp[10];
   return(&temp[11]); // returns a pointer to the local int
}

void Victim() {
   int* ptr;
   ptr = TAB();
   *ptr = 42;
    fprintf(stderr,"ptr =%d\n", *ptr);
}

int main(int argc, char *argv[]) {
    Victim();
    return 0;
}

