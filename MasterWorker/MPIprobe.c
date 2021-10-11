#include <assert.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

int mpi_handleModeling(float *data, int size_s, float *modResult, int size_i);

int main(int argc, char *argv[])
{
	int     pe, root_pe=0, npes, end_of_work, one_shot;
	int     nt; 
	int     flag, itrace, fldr,kpe, work, i;
	int     shot_tag, size_i, size_s;
	float   scl, modResult[150], data[100];
	MPI_Status status;
	MPI_Request request;

	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &npes );
	MPI_Comm_rank( MPI_COMM_WORLD, &pe );

	if (npes <= 1) {
		fprintf(stderr, "Number of processors is smaller than 2, use a minimum of 2 processors !\n");
		MPI_Finalize();
		exit(0);
	}
	
/* =========== Open input parameter file  ========== */

	if (pe == root_pe) {

	/* ..... parameter file read by master only */
		size_i=150;
		scl=5.0;
	
	}
	size_s=100;				

	/* communicate  information to other PE's */

    MPI_Bcast( &size_i, 1, MPI_INT, root_pe, MPI_COMM_WORLD );
    MPI_Bcast( &scl, 1, MPI_FLOAT, root_pe, MPI_COMM_WORLD );

/* =============== forward modeling over all shot positions ================= */

	kpe=1;
	end_of_work = 0;
	while (!end_of_work) {

		work = mpi_handleModeling(data, size_s, modResult, size_i);

		if (work) { /* only the working pe's will do this */
			/* start forward modeling of data */
			for (i=0;i<size_i; i++) modResult[i]=(float)work;

			sleep(1);
			fprintf(stderr, "*** Source position %d processed by pe %d ***\n", kpe++, pe);

		} /* end of processing */
		else {
			end_of_work = 1;
		}

	} /* end of data while loop */


	fprintf(stderr, "pe %d: before last barrier\n", pe);
	MPI_Barrier(MPI_COMM_WORLD);
	
	/* =============== continue with scheme  ================= */

	
	MPI_Finalize();

	return 0;
}




int mpi_handleModeling(float *data, int size_s, float *modResult, int size_i)
{       
    static int first=1, work_done=0;
    int workingpes, root_pe=0, pe, npes, flag, i, source, work;
	float bufResult[150];
	
	static float *bufex;
	int data_work_tag, data_request_tag, data_tag, data_xrcv_tag, data_hdr_tag;
	int src_nx_tag, src_request_tag, src_hdr_tag, src_tag, image_tag, shot_tag;
    MPI_Status status;
    static MPI_Request reqDd, reqIm; 
                
    MPI_Comm_size( MPI_COMM_WORLD, &npes );
    MPI_Comm_rank( MPI_COMM_WORLD, &pe );

	data_request_tag = 1;
	data_work_tag    = 2;
	data_tag         = 3;
	data_xrcv_tag    = 4;
	data_hdr_tag     = 5;
	src_request_tag  = 6;
	src_nx_tag       = 7;
	src_tag          = 8;
	image_tag        = 9;
	shot_tag         = 10;

	work=50; /* just for testing there are 50 tasks to do */
	for (i=0;i<size_s; i++) data[i]=(float)work;
	if (pe == root_pe) {
		if (first) { /* if needed do some initialisations only needed for the first time */
			first = 0;
		}
		workingpes=npes-1;
		while (workingpes) {

			/* root_pe waiting for request for data */
			fprintf(stderr,"*0* root_pe waiting for shot data request\n");
			MPI_Probe(MPI_ANY_SOURCE, data_request_tag, MPI_COMM_WORLD, &status);
			source = status.MPI_SOURCE;
			MPI_Recv(&flag, 1, MPI_INT, source, status.MPI_TAG, MPI_COMM_WORLD, &status);
			fprintf(stderr,"*0* root_pe received request from %d for sending work=%d \n", source, work);
	
			/* check if there is more work todo  */
			if (work != 0) { /* there is more work send it to processor source */

				MPI_Send(&work, 1, MPI_INT, source, data_work_tag, MPI_COMM_WORLD);
				
				MPI_Isend(data, size_s, MPI_FLOAT, source, data_tag, MPI_COMM_WORLD, &reqDd);
			}
			else { /* no more work to do this processor can stop */
				work = 0;
				MPI_Send(&work, 1, MPI_INT, source, data_work_tag, MPI_COMM_WORLD);
				workingpes--;
				fprintf(stderr,"pe %d can leave the loop, working pe's left %d\n", source, workingpes);
			}

			/* receive modeling results from working PE's */
			
			MPI_Iprobe(MPI_ANY_SOURCE, shot_tag, MPI_COMM_WORLD, &flag, &status);
			while(flag) { /* receive result from source */
				size_i=150;
				MPI_Recv(&modResult[0], size_i, MPI_FLOAT, MPI_ANY_SOURCE, shot_tag, MPI_COMM_WORLD, &status);
				fprintf(stderr,"master received the %f results from pe %d \n", modResult[0], status.MPI_SOURCE);

				/* process the results of this forward modeled data */
				
				MPI_Iprobe(MPI_ANY_SOURCE, shot_tag, MPI_COMM_WORLD, &flag, &status);
			}

			/* define work task */
			if (work) {
				MPI_Wait(&reqDd, MPI_STATUS_IGNORE); /* make sure previous task has been sent */
				work--; /* one more task less to do */
				for (i=0;i<size_s; i++) data[i]=(float)work;
			}

		} /* end of workingpes loop */
		
		fprintf(stderr,"Master left the workingpes loop\n");


	} /* end of root part */
	else {  /* non root part */
		/* sent calculated model back to root */
		/* blocking send to root_pe */
		if (work_done) {
			size_i = 150;
			memcpy(&bufResult[0], &modResult[0], size_i*sizeof(float));

			fprintf(stderr,"#1# pe %d sending modeled data of source %f to root\n", pe, modResult[0]);
			MPI_Send(&bufResult[0], size_i, MPI_FLOAT, root_pe, shot_tag, MPI_COMM_WORLD);
		}

		/* request new work from root pe */
		work = 0;
		MPI_Send(&work, 1, MPI_INT, root_pe, data_request_tag, MPI_COMM_WORLD);
		MPI_Recv(&work, 1, MPI_INT, root_pe, data_work_tag, MPI_COMM_WORLD, &status);
		
		/* there is no more work to be done: leave the while loop */
		if (work == 0) { 
			fprintf(stderr,"#0# pe %d has no more work: leaving while loop\n", pe);
			return 0;
		}
		else { /* there is more work to do */
			size_s = 100;
			MPI_Recv(data, size_s, MPI_FLOAT, root_pe, data_tag, MPI_COMM_WORLD, &status);
			fprintf(stderr,"#0# pe %d has received more work for %f \n", pe, data[0]);

			sleep(1);

			work_done++;
		}

	}


	return work;
}

