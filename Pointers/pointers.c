#include<stdlib.h>
#include<stdio.h>

int main (int argc, char *argv[])
{
	unsigned int *ptr, p, q, i;
	size_t size;
	float *A;

	size = 1024*1024*1024;
//	size = 1024*1024;
	ptr = (unsigned int *)malloc(size*sizeof(unsigned int));
	ptr[0] = 123;
	p = 321;
	A = (float *)malloc(0100*sizeof(float));
/* valid assignment only on 32 address-space */
	q = ptr;
	fprintf(stderr,"sizeof(size_t int)=%d\n", sizeof(size_t));
	fprintf(stderr,"sizeof(unsigned int)=%d\n", sizeof(unsigned int));
	fprintf(stderr,"sizeof(ptr)=%d\n", sizeof(ptr));
	fprintf(stderr,"sizeof(*ptr)=%d\n", sizeof(*ptr));
	fprintf(stderr,"ptr=%d ptr[0]=%d *ptr=%d p=%d q=%d\n",ptr, ptr[0], *ptr, p,  q);
	fprintf(stderr,"ptr=0x%016llX ptr[0]=0x%X *ptr=0x%X p=0x%X q=0x%016llX\n",ptr, ptr[0], *ptr, p,  q);
	fprintf(stderr,"ptr=0x%016lX ptr[0]=0x%X *ptr=0x%X p=0x%X q=0x%016lX\n",ptr, ptr[0], *ptr, p,  q);
	for (i=0; i<100;i++) A[i] = (float)i;
	return 0;
}
