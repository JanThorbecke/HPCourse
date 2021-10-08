#include <fcntl.h>
#include <unistd.h>
#include <signal.h>
#include <inttypes.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
/*#include <aio.h>*/
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>

char *filename;
void write_data_unbuf(float *data, int size, int chuncksize);
void read_data_unbuf(float *data, int size, int chuncksize);

void write_data_buf(float *data, int size, int chuncksize);
void read_data_buf(float *data, int size, int chuncksize);

void write_data_mmap(float *data, int size);
void read_data_mmap(float *data, int size);

void write_data_ascii(float *data, int size, int chuncksize);
void read_data_ascii(float *data, int size, int chuncksize);

void write_data_aio(float *data, int size);
void read_data_aio(float *data, int size);
static void sig_usr1(int signo);

void putgarbage(int size);
float wallclock_time();

/*********************************
* Examples of IO:
*   - buffered (fopen ...)
*   - unbuffered (open ...)
*   - a-synchronous
*   - memory mapped
************************************/

int main(int argc, char **argv)
{
	int fd, retval, i;
	size_t size, chuncksize;
	int ch;
	float t0, t1, t2, t3;
	float *data;
	struct stat status;
	char syscall[100];

	size = 4*1024*1024; // size in number of floats in array data
  /* cygwin users change size in line above to 32*1024*1024 */
	filename = "/tmp/fileio.bin";
	while ((ch = getopt(argc, argv, "b:f:")) != -1) {
		switch (ch) {
			case 'b':
				size = 1024*atoi(optarg)/sizeof(float);
				break;
			case 'f':
				filename = optarg;
				break;
			default:
				break;
		}
	}
	fprintf(stderr,"filename=%s size=%ld KB\n", filename, size);
	if (size < 2048) {
		fprintf(stderr,"The outpute file size (-b option) must be larger or equal to 8\n");
		exit(0);
	}

	if ( (fd = open(filename, 
		(O_RDWR | O_TRUNC | O_CREAT), 0644 )) < 0){
		perror("main"); 
		exit(1);
	}
	if((retval = fstat(fd, &status) != 0)) {
		perror("fstat failed");
		exit(1);
	}
	if (close(fd) < 0) {
		perror("main");
		exit(1);
	}

	fprintf(stderr,"prefered block size for IO    = %d bytes\n", 
		status.st_blksize);
	fprintf(stderr,"system pagesize               = %d bytes\n", 
		getpagesize() );
	fprintf(stderr,"Data size = %.3f MB\n", size/(1024.0*1024.0));

/* generate data to write to file */

	chuncksize = 8192;
	data = (float *)malloc((2*size)*sizeof(float));
	for (i = 0; i < size; i++) data[i] = (float) i;
	size *= sizeof(float); // size in Bytes

/* Write data to file */

	fprintf(stderr,"***** WRITING DATA *****\n");


	fprintf(stderr,"Kernel open/write IO, with sync\n");
	while (chuncksize <= size) {
		write_data_unbuf(data, size, chuncksize);
		chuncksize *= 2;
	}
	sleep(3);

	fprintf(stderr,"ANSI C fopen/fwrite IO\n");
	chuncksize = 8192;
	while (chuncksize <= size) {
		write_data_buf(data, size, chuncksize);
		chuncksize *= 2;
	}
	sleep(3);

	fprintf(stderr,"Memory mapped IO\n");
	write_data_mmap(data, size);


/* put garbage in kernel filesystem buffer */

	putgarbage(4096*128);
	fprintf(stderr,"***** READING DATA *****\n");

/* Read data from file */

	fprintf(stderr,"Kernel open/read IO, with sync\n");
	chuncksize = 8192;
	while (chuncksize <= size) {
		read_data_unbuf(data, size, chuncksize);
		chuncksize *= 2;
		putgarbage(4096*128);
	}
	sleep(3);

	fprintf(stderr,"ANSI C fopen/fread IO\n");
	chuncksize = 8192;
	while (chuncksize <= size) {
		read_data_buf(data, size, chuncksize);
		chuncksize *= 2;
		putgarbage(4096*128);
	}
	sleep(3);

	fprintf(stderr,"Memory mapped IO\n");
	read_data_mmap(data, size);

	
	fprintf(stderr,"*****  ASCII DATA *****\n");

	chuncksize = size;
	while (chuncksize <= size) {
		write_data_ascii(data, size, chuncksize);
		chuncksize *= 2;
	}
	sleep(3);
	
	chuncksize = size;
	while (chuncksize <= size) {
		read_data_ascii(data, size, chuncksize);
		chuncksize *= 2;
	}
	sleep(3);
	
	/* remove data file */
	sprintf(syscall,"rm %s\n", filename);
	system(syscall);
	exit(0);
}

/***************************************************
* Use system kernels open, read and write
***************************************************/

void write_data_unbuf(float *data, int size, int chuncksize)
{
	int fd, n, N, esz;
	float t0, t1, t;
	void *p;

	if ( (fd = open(filename, (O_RDWR | O_TRUNC | O_CREAT), 
		(S_IRUSR | S_IWUSR))) < 0){
		perror("write_data_unbuf"); 
		exit(1);
	}

	N = 0;
	p = (void *)data;
	esz = sizeof(*data);
	t0 = wallclock_time();
	while ( (n = write(fd, p, chuncksize)) > 0 && N < size) {
		N += n;
		p = data+N/esz;
	}
	sync();
	t1 = wallclock_time();
	t = t1 - t0;

	fprintf(stderr,"written %6.3f MB : %7d loops : chuncksize %7d : %.3f s. => %.3f MB/s\n", N/1e6, N/chuncksize, chuncksize, t, size/(1e6*t));
		
	if (close(fd) < 0) {
		perror("write_data_unbuf");
		exit(1);
	}

	return;
}

void read_data_unbuf(float *data, int size, int chuncksize)
{
	int fd, n, N, esz;
	float t0, t1, t;
	void *p;

	if ( (fd = open(filename, O_RDWR, (S_IRUSR | S_IWUSR))) < 0){
		perror("read_data_unbuf"); 
		exit(1);
	}

	N = 0;
	p = (void *)data;
	esz = sizeof(*data);
	t0 = wallclock_time();
	while ( (n = read(fd, p, chuncksize)) > 0 && N < size) {
		N += n;
		p = data+N/esz;
	}
	sync();
	t1 = wallclock_time();
	t = t1 - t0;

	fprintf(stderr,"read %6.3f MB : %7d loops : chuncksize %7d : %.3f s. => %.3f MB/s\n", N/1e6, N/chuncksize, chuncksize, t, size/(1e6*t));
		
	if (close(fd) < 0) {
		perror("read_data_unbuf");
		exit(1);
	}

	return;
}

/**************************************************/

void write_data_buf(float *data, int size, int chuncksize)
{
	FILE *fp;
	int n, N, esz;
	float t0, t1, t;
	char *buffer;
	void *p;

	if ( (fp = fopen(filename, "w")) == NULL) {
		perror("write_data_buf"); 
		exit(1);
	}

/* _IONBF := No Buffer; _IOFBF := Full Buffer */

/*
	buffer = (char *)calloc(chuncksize, sizeof(float));
	setvbuf(fp, buffer, _IONBF, chuncksize*sizeof(float));
*/
//	setbuf(fp, NULL);

	N = 0;
	p = (void *)data;
	esz = sizeof(*data);
	t0 = wallclock_time();
	while ( (n = fwrite(p, 1, chuncksize, fp)) > 0 && N < size) {
		N += n;
		p = data+N/esz;
	}
	t1 = wallclock_time();
	t = t1 - t0;

	fprintf(stderr,"written %6.3f MB : %7d loops : chuncksize %7d : %.3f s. => %.3f MB/s\n", N/1e6, N/chuncksize, chuncksize, t, size/(1e6*t));
		
	if (fclose(fp) < 0) {
		perror("write_data_buf");
		exit(1);
	}

/*
	free(buffer);
*/
	return;
}

void read_data_buf(float *data, int size, int chuncksize)
{
	FILE *fp;
	int n, N, esz;
	float t0, t1, t;
	char *buffer;
	void *p;

	if ( (fp = fopen(filename, "r")) == NULL) {
		perror("read_data_buf"); 
		exit(1);
	}

//	setbuf(fp, NULL);

	N = 0;
	p = (void *)data;
	esz = sizeof(*data);
	t0 = wallclock_time();
	while ( (n = fread(p, 1, chuncksize, fp)) > 0 && N < size) {
		N += n;
		p = data+N/esz;
	}
	t1 = wallclock_time();
	t = t1 - t0;

	fprintf(stderr,"read %6.3f MB : %7d loops : chuncksize %7d : %.3f s. => %.3f MB/s\n", N/1e6, N/chuncksize, chuncksize, t, size/(1e6*t));
		
	if (fclose(fp) < 0) {
		perror("read_data_buf");
		exit(1);
	}

/*
	free(buffer);
*/
	return;
}
/*******************************************
* memory mapped IO
*******************************************/

void write_data_mmap(float *data, int size)
{
	int fd;
	float t0, t1, t;
	void *address;

	if ( (fd = open(filename, O_RDWR, (S_IRUSR | S_IWUSR))) < 0){
		perror("write_data_mmap"); 
		exit(1);
	}
	lseek(fd,0,SEEK_SET);

	t0 = wallclock_time();
	address = mmap(NULL, size, (PROT_READ | PROT_WRITE), MAP_PRIVATE, fd, 0);
	if (address == MAP_FAILED) perror("write_data_mmap");
	memcpy(address, data, size);
//	sync();
	t1 = wallclock_time();
	t = t1 - t0;

	fprintf(stderr,"written %6.3f MB : %7d loops : chuncksize %7d : %.3f s. => %.3f MB/s\n", size/1e6, size/size, size, t, size/(1e6*t));
		
	if (close(fd) < 0) {
		perror("write_data_mmap");
		exit(1);
	}

	return;
}

void read_data_mmap(float *data, int size)
{
	int fd;
	float t0, t1, t;
	void *address;

	if ( (fd = open(filename, O_RDWR, (S_IRUSR | S_IWUSR))) < 0){
		perror("read_data_mmap"); 
		exit(1);
	}
	lseek(fd,0,SEEK_SET);

	t0 = wallclock_time();
	address = mmap(NULL, size, (PROT_READ | PROT_WRITE), MAP_PRIVATE, fd, 0);
	if (address == MAP_FAILED) perror("write_data_mmap");
	memcpy(data, address, size);
//	sync();
	t1 = wallclock_time();
	t = t1 - t0;

	fprintf(stderr,"read %6.3f MB : %7d loops : chuncksize %7d : %.3f s. => %.3f MB/s\n", size/1e6, size/size, size, t, size/(1e6*t));
		
	if (close(fd) < 0) {
		perror("read_data_mmap");
		exit(1);
	}

	return;
}

/*******************************************
 * ASCII
 *******************************************/

void write_data_ascii(float *data, int size, int chuncksize)
{
	FILE *fp;
	int n, N, esz;
	float t0, t1, t;
	char *buffer;
	void *p;
	
	if ( (fp = fopen(filename, "w")) == NULL) {
		perror("write_data_ascii"); 
		exit(1);
	}
		
	N = 0;
	esz = sizeof(*data);
	t0 = wallclock_time();
	while ( N < size ) {
		fprintf(fp,"%f\n",data[N/esz]);
		N += esz;
	}
	fflush(fp);
	t1 = wallclock_time();
	t = t1 - t0;
	
	fprintf(stderr,"written %6.3f MB : %7d loops : chuncksize %7d : %.3f s. => %.3f MB/s\n", N/1e6, N/chuncksize, chuncksize, t, size/(1e6*t));
	
	if (fclose(fp) < 0) {
		perror("write_data_ascii");
		exit(1);
	}
	return;
}


void read_data_ascii(float *data, int size, int chuncksize)
{
	FILE *fp;
	int n, N, esz;
	float t0, t1, t;
	char *buffer;
	void *p;
	   
	if ( (fp = fopen(filename, "r")) == NULL) {
	   perror("read_data_buf"); 
	   exit(1);
	}
	   
	N = 0;
	esz = sizeof(float);
	t0 = wallclock_time();
	while ( N < size ) {
		   fscanf(fp,"%f\n",&data[N/esz]);
		   N += esz;
	}
	fflush(fp);
	t1 = wallclock_time();
				  
	t = t1 - t0;
	   
	fprintf(stderr,"read %6.3f MB : %7d loops : chuncksize %7d : %.3f s. => %.3f MB/s\n", N/1e6, N/chuncksize, chuncksize, t, size/(1e6*t));
	   
	if (fclose(fp) < 0) {
		perror("read_data_ascii");
		exit(1);
	}
	   
	return;
}
		   
		   
		   


/*******************************************
* asynchronous IOO
*******************************************/
/*
void write_data_aio(float *data, int size)
{
	aioinit_t init;
	aiocb_t aiocbp;
	int fd;
	float t0, t1, t;
	void *address;

	init.aio_threads = 5;
	init.aio_locks = 3;
	init.aio_numusers = 5;
	aio_sgi_init(&init);

	if ( (fd = open(filename, (O_RDWR | O_TRUNC | O_CREAT ), 
		0644 )) < 0){
		perror("write_data_mmap"); 
		exit(1);
	}

	aiocbp.aio_nbytes = size;
	aiocbp.aio_fildes = fd;
	aiocbp.aio_buf = (float *)data;	
	aiocbp.aio_sigevent.sigev_notify = SIGEV_SIGNAL;
	aiocbp.aio_sigevent.sigev_signo = SIGUSR1;
	
	if (signal(SIGUSR1, sig_usr1) == SIG_ERR)
		fprintf(stderr,"signal(SIGUSR1) error\n");


	t0 = wallclock_time();

	if (aio_write(&aiocbp) < 0)
	aio_return(&aiocbp);
	sleep(10);

	t1 = wallclock_time();
	t = t1 - t0;

	fprintf(stderr,"written %6.3f MB : %.3f s. => %.3f MB/s\n", size/1e6, t, size/(1e6*t));
		
	return;
}

static void sig_usr1(int signo)
{
	fprintf(stderr,"SIGUSR1 signal\n");

	return;
}

void read_data_aio(float *data, int size)
{
	int fd;
	float t0, t1, t;
	void *address;

	if ( (fd = open(filename, O_RDWR)) < 0){
		perror("read_data_mmap"); 
		exit(1);
	}
	lseek(fd,0,SEEK_SET);

	t0 = wallclock_time();
	address = mmap(NULL, size, (PROT_READ | PROT_WRITE), MAP_PRIVATE, fd, 0);
	if (address == MAP_FAILED) perror("write_data_mmap");
	memcpy(data, address, size);
	sync();
	t1 = wallclock_time();
	t = t1 - t0;

	fprintf(stderr,"read %6.3f MB : %7d loops : chuncksize %7d : %.3f s. => %.3f MB/s\n", size/1e6, size/size, size, t, size/(1e6*t));
		
	if (close(fd) < 0) {
		perror("read_data_mmap");
		exit(1);
	}

	return;
}
*/
/*********************************
* put garbage in filesytem buffer 
*********************************/

void putgarbage(int size)
{
	int fd, n;
	float *nep;

	if ( (fd = open("/tmp/garbage", 
		(O_RDWR | O_TRUNC | O_CREAT), 0644 )) < 0){
		perror("main"); 
		exit(1);
	}
	lseek(fd,100*size,SEEK_SET);

	nep = (float *)calloc(size, sizeof(float));
	n = write(fd, nep, size);
	assert(n == size);
	if (close(fd) < 0) {
		perror("main");
		exit(1);
	}
	system("rm /tmp/garbage\n");
	return;
}

/*
* function to calculate wallclock-time
*/

float wallclock_time()
{
	struct timeval s_val;
	static struct timeval b_val;
	float time;
	static int base=0;

	gettimeofday(&s_val,0);

	if (!base) {
		b_val = s_val;
		base = 1;
		return 0.0;
	}

	time = (float)(s_val.tv_sec-b_val.tv_sec) + 
		   (float)(1e-6*(s_val.tv_usec-b_val.tv_usec));

	return (float)time;
}

