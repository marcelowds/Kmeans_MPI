#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
#define DIM 3

int main(int argc, char **argv) {
	int my_rank,n_procs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
	int kn[2];
	int i,c,k,j, n,*cluster,*count,*global_count,flips,global_flips,numit,color;
	double *x, *x_i,*mean,*sum,*global_sum,dx,dmin;
	int ret;
	if(my_rank==0){
		ret=scanf("%d", &kn[0]);
		ret=scanf("%d", &kn[1]);
	}
	MPI_Bcast(kn, 2, MPI_INT, 0, MPI_COMM_WORLD);
	k=kn[0];n=kn[1];
	int tam_p_proc=n/n_procs+1;
	x = (double *)malloc(sizeof(double)*DIM*tam_p_proc*n_procs);

	int min, max;
	min=tam_p_proc*my_rank;
	if(my_rank!=n_procs-1){max=tam_p_proc*(my_rank+1);}else{max=n;}

	x_i = (double *)malloc(sizeof(double)*DIM*(tam_p_proc));
	mean = (double *)malloc(sizeof(double)*DIM*k);
	sum= (double *)malloc(sizeof(double)*DIM*k);
	global_sum= (double *)malloc(sizeof(double)*DIM*k);
	cluster = (int *)malloc(sizeof(int)*n);
	count = (int *)malloc(sizeof(int)*k);
	global_count = (int *)malloc(sizeof(int)*k);

	if(my_rank==0){
		for (i = 0; i<k; i++)
			ret=scanf("%lf %lf %lf", mean+i*DIM, mean+i*DIM+1, mean+i*DIM+2);
		for (i = 0; i<n; i++)
			ret=scanf("%lf %lf %lf", x+i*DIM, x+i*DIM+1, x+i*DIM+2);
	}

	MPI_Bcast(mean, DIM*k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(x, DIM*tam_p_proc, MPI_DOUBLE, x_i,DIM*tam_p_proc,MPI_DOUBLE, 0, MPI_COMM_WORLD);
	for (i = 0; i<n; i++){cluster[i] = 0;}
	free(x);
	global_flips=n;
	numit=0;

	while (global_flips>0) {
		numit++;
		flips = 0;
		global_flips=0;
		for (j = 0; j < k; j++) {
			count[j] = 0;
			for (i = 0; i < DIM; i++){
				sum[j*DIM+i] = 0.0;
				global_sum[j*DIM+i] = 0.0;
			}
		}
		for (i = min; i < max; i++) {
			dmin = -1; color = cluster[i];
			for (c = 0; c < k; c++) {
				dx = 0.0;
				for (j = 0; j < DIM; j++)
					dx +=  (x_i[(i-min)*DIM+j] - mean[c*DIM+j])*(x_i[(i-min)*DIM+j] - mean[c*DIM+j]);
				if (dx < dmin || dmin == -1) {
					color = c;
					dmin = dx;
				}
			}
			if (cluster[i] != color) {
				flips++;
				cluster[i] = color;
		  	}

		}
		for (i = min; i < max; i++) {
			count[cluster[i]]++;
			for (j = 0; j < DIM; j++){
				sum[cluster[i]*DIM+j] += x_i[(i-min)*DIM+j];
			}

		}
		MPI_Allreduce(count, global_count,k,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(sum, global_sum,k*DIM,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(&flips, &global_flips,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		for (i = 0; i < k; i++) {
			for (j = 0; j < DIM; j++) {
				mean[i*DIM+j] = global_sum[i*DIM+j]/global_count[i];
  			}
		}
	}
	ret=ret+0;
	if(my_rank==0){
		for (i = 0; i < k; i++) {
			for (j = 0; j < DIM; j++)
				printf("%5.2f ", mean[i*DIM+j]);
			printf("\n");
		}
	}
	MPI_Finalize();
	return(0);
}
