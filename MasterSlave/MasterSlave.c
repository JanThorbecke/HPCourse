#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#define WORKTAG 1
#define DIETAG 2

/* Local functions */

#define unit_of_work_t int
static void master(void);
static void slave(void);
static unit_of_work_t get_next_work_item(void);
static void process_results(double result);
static double do_work(unit_of_work_t work);


int
main(int argc, char **argv)
{
  int myrank;

  /* Initialize MPI */

  MPI_Init(&argc, &argv);

  /* Find out my identity in the default communicator */

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0) {
    master();
  } else {
    slave();
  }

  /* Shut down MPI */

  MPI_Finalize();
  return 0;
}


static void
master(void)
{
  int ntasks, rank;
  unit_of_work_t work;
  double result;
  MPI_Status status;

  /* Find out how many processes there are in the default
     communicator */

  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  /* Seed the slaves; send one unit of work to each slave. */

  for (rank = 1; rank < ntasks; ++rank) {

    /* Find the next item of work to do */

    work = get_next_work_item();

    /* Send it to each rank */

	fprintf(stderr,"Master sending work %d to worker %d\n", work, rank);
    MPI_Send(&work,             /* message buffer */
             1,                 /* one data item */
             MPI_INT,           /* data item is an integer */
             rank,              /* destination process rank */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */
  }

  /* Loop over getting new work requests until there is no more work
     to be done */

  work = get_next_work_item();

  while (work != 0) {

    /* Receive results from a slave */

    MPI_Recv(&result,           /* message buffer */
             1,                 /* one data item */
             MPI_DOUBLE,        /* of type double real */
             MPI_ANY_SOURCE,    /* receive from any sender */
             MPI_ANY_TAG,       /* any type of message */
             MPI_COMM_WORLD,    /* default communicator */
             &status);          /* info about the received message */

    /* Send the slave a new work unit */
	fprintf(stderr,"Master received result from worker %d\n", status.MPI_SOURCE);
	fprintf(stderr,"Master sending work %d to worker %d\n", work, status.MPI_SOURCE);

    MPI_Send(&work,             /* message buffer */
             1,                 /* one data item */
             MPI_INT,           /* data item is an integer */
             status.MPI_SOURCE, /* to who we just received from */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */

    /* Get the next unit of work to be done */

    work = get_next_work_item();
  }

  /* There's no more work to be done, so receive all the outstanding
     results from the slaves. */

  for (rank = 1; rank < ntasks; ++rank) {
    MPI_Recv(&result, 1, MPI_DOUBLE, MPI_ANY_SOURCE,
             MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }

  /* Tell all the slaves to exit by sending an empty message with the
     DIETAG. */

  for (rank = 1; rank < ntasks; ++rank) {
    MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
  }
}


static void 
slave(void)
{
  unit_of_work_t work;
  double result;
  MPI_Status status;
  int myrank;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  while (1) {

    /* Receive a message from the master */

    MPI_Recv(&work, 1, MPI_INT, 0, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);

	fprintf(stderr,"Worker %d received from master task %d\n", myrank, work);
    /* Check the tag of the received message. */

    if (status.MPI_TAG == DIETAG) {
		fprintf(stderr,"Worker %d received no more work from master\n", myrank);
      return;
    }

    /* Do the work */

    result = (double)do_work(work);

    /* Send the result back */

    MPI_Send(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
}


static unit_of_work_t 
get_next_work_item(void)
{
	static int tobedone=100;
  /* Fill in with whatever is relevant to obtain a new unit of work
     suitable to be given to a slave. */
	tobedone--;
	return tobedone;
}


static void 
process_results(double result)
{
  /* Fill in with whatever is relevant to process the results returned
     by the slave */
}


static double
do_work(unit_of_work_t work)
{
  /* Fill in with whatever is necessary to process the work and
     generate a result */
	sleep(1);
	return 1;
}

