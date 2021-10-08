#include<stdlib.h>
#include<stdio.h>

int main (int argc, char *argv[])
{
	int arr[10], j, *ptr;

	arr[0] = 7;
	ptr = &arr[0];
	*(ptr+2)=10;
	fprintf(stderr,"arr[0]=%d arr[1]=%d arr[2]=%d prt=%x\n", arr[0], arr[1], arr[2], arr);

	return 0;
}
