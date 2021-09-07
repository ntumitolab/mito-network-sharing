// emergence of effective genomes on mitochondrial encounter networks

// analysed as the "bingo" game where nodes must collect different tokens by exchanging with their neighbours on a network

// max and mean bingo statistics for a given adjacency matrix, and for negative controls based on those graph stats
// there are lots of these controls. ER networks, SF networks constructed using various protocols, "cliquey" networks likewise, and networks arising from basic physical simulation
// compile with igraph library e.g. gcc bingo.c -I /usr/local/include/igraph -ligraph -lm -g -o bingo.ce
// OR, to avoid using igraph, comment out the line #define _USEIGRAPH_ 1

#define _USEIGRAPH_ 1

#ifdef _USEIGRAPH_
#include <igraph.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int OUTPUT_SIMULATION_TRACES = 0;  // output simulated mitochondrial trajectories?
int REMOVE_SINGLETONS = 0;         // remove singletons from bio networks?
int OUTPUT_TRAJECTORIES = 0;       // number of sims from which to store time-series output 
int INACTIVITY = 0;                // whether simulated mitos can be inactive

#define BINGOTYPE 0                // if there is time ordering of edges (as in bio and simulation), 0 uses early edges first, 1 uses random edge ordering
#define MAXN 10000
#define MAXE 10000

#define NIT 10         // number of random iterations to run (was 100)
#define SAMPLING 1000
#define SAMPLINGIT 10

// for Mac and Linux
#define RND drand48()
// for Windows
// #define RND ((double) rand() / (RAND_MAX))

// produce gaussian random number
double gsl_ran_gaussian(const double sigma)
{
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = -1 + 2 * RND;
      y = -1 + 2 * RND;

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
}

// read in adjacency matrix from file
// no error checking or anything -- not robust
void AMFromFile(char *fstr, int *amlist, int *amn, int *nnodes)
{
  FILE *fp;
  int tmp1, tmp2, tmp3;
  char ch;
  int i, j, k;
  int tmp;
  int degree;
  int singles[MAXN];
  int nsingles;
  int min, max;
  
  fp = fopen(fstr, "r");
  if(fp == NULL)
    {
      printf("Couldn't open %s\n", fstr);
      exit(0);
    }
  *amn = 0; max = 0;
  // skip header line
  do{ ch = fgetc(fp); }while(ch != '\n');
  do{
    // skip first entry in each row (frame detail), record second and third (node labels)
    fscanf(fp, "%i,%i,%i", &tmp1, &tmp2, &tmp3);
    if(feof(fp)) break;
    amlist[(*amn)++] = tmp2;
    amlist[(*amn)++] = tmp3;

    // am list file ends with a fake self-edge describing the total number of nodes in the network
    if(tmp2 == tmp3)
      {
	max = tmp2;
	(*amn) -= 2;
      }
    
  }while(!feof(fp));

  min = 9999;
  for(j = 0; j < *amn; j++)
    {
      if(amlist[j] < min)
	min = amlist[j];
    }

  //  printf("  Prior to removing singletons and normalising, minimum node is %i and maximum node is %i\n", min, max);

  if(REMOVE_SINGLETONS)
    {
      nsingles = 0;
      for(i = min; i <= max; i++)
	{
	  degree = 0;
	  for(j = 0; j < *amn; j++)
	    degree += (amlist[j] == i);
	  if(degree == 0)
	    singles[nsingles++] = i;
	}
      //      printf("  %i singletons\n", nsingles);
      
      for(j = 0; j < *amn; j++)
	{
	  tmp = amlist[j];
	  for(k = 0; k < nsingles; k++)
	    {
	      if(singles[k] <= amlist[j])
		tmp--;
	    }
	  amlist[j] = tmp;
	}
      (max) = 0;
      for(j = 0; j < *amn; j++)
	{
	  if(amlist[j] > max)
	    max = amlist[j];
	}
    }

  min = 9999;
  for(j = 0; j < *amn; j++)
    {
      if(amlist[j] < min)
	min = amlist[j];
    }
  
  for(j = 0; j < *amn; j++)
    amlist[j] -= min;
  (max) -= min;
  min -= min;
  
  /*  if(REMOVE_SINGLETONS)
      printf("  I read minimum node number %i and maximum node number %i\n", min, max);
      else
      printf("  I read minimum node number %i and top node label %i\n", min, max);*/

  (*nnodes) = max+1;
  (*amn) /= 2;
  fclose(fp);
}

// get minimum and maximum integers from a file
// assumes no negative integers
/*void StatsFromFile(char *fstr, int *min, int *max)
  {
  FILE *fp;
  int tmp1, tmp2, tmp3;
  char ch;
  
  // this has been changed to suit the new I/O format. node labels start from 1 and run to n
  // the final entry in the am list file is a fake self-edge on the nth node
  
  *min = 1;

  fp = fopen(fstr, "r");
  if(fp == NULL)
  {
  printf("Couldn't open %s\n", fstr);
  exit(0);
  }
  // skip header line
  do{ ch = fgetc(fp); }while(ch != '\n');
  do{
  fscanf(fp, "%i,%i,%i", &tmp1, &tmp2, &tmp3);
  if(feof(fp)) break;
  *max = tmp3;
  }while(!feof(fp));
  fclose(fp);
  }*/

// build Erdos-Renyi network matching n and e statistics
void AMErdosRenyi(int n, int e, int *amlist, int *amn)
{
  int j;
  int r1, r2;
  int *am;

  am = (int*)malloc(sizeof(int)*n*n);
  for(j = 0; j < n*n; j++)
    am[j] = 0;
  
  for(j = 0; j < e; j++)
    {
      r1 = RND*n;
      do{ r2 = RND*n; }while(r1 == r2 && am[r1*n+r2] == 0);
      amlist[2*j+0] = r1;
      amlist[2*j+1] = r2;
      am[r1*n+r2] = am[r2*n+r1] = 1;
    }
  *amn = e;
  free(am);
}

// build cliquey graph with cliques of size c
// flavour == 0, disconnected cliques
// flavour == 1, single connection between cliques
// flavour == 2, pad edges with random connections after building cliques
// flavour == 3, flavours 1+2
void AMCliquey(int n, int e, int *amlist, int *amn, int c, int flavour)
{
  int i, j, k;
  int r1, r2;
  int nclique_edges, nclique_nodes, nclique;

  nclique_edges = floor(e / (c*(c-1)/2 + 1));
  nclique_nodes = floor(n / c);
  if(nclique_edges < nclique_nodes) nclique = nclique_edges;
  else nclique = nclique_nodes;
  
  (*amn) = 0;
  for(i = 0; i < c*nclique; i += c)
    {
      for(j = i; j < i+c; j++)
	{
	  for(k = j+1; k < i+c; k++)
	    {
	      amlist[2*(*amn)+0] = j;
	      amlist[2*(*amn)+1] = k;
	      (*amn)++;
	    }
	}
    }
  if(flavour == 1 || flavour == 3)
    {
      for(i = c; i < c*nclique; i += c)
	{
	  amlist[2*(*amn)+0] = i;
	  amlist[2*(*amn)+1] = i-c;
	  (*amn)++;
	}
    }
  if(flavour == 2 || flavour == 3)
    {
      for(j = *amn; j < e-1; j++)
	{
	  amlist[2*j+0] = RND*n;
	  amlist[2*j+1] = RND*n;
	}
      *amn = j;
    }
  if((*amn) >= e) *amn = e;
}

// build star graph
void AMStar(int n, int e, int *amlist, int *amn)
{
  int i;
  
  (*amn) = 0;
  for(i = 1; i < n; i++)
    {
      amlist[2*(*amn)+0] = 0;
      amlist[2*(*amn)+1] = i;
      (*amn)++;
    }      
  for(i = *amn; i < e-1; i++)
    {
      amlist[2*i+0] = RND*n;
      amlist[2*i+1] = RND*n;
    }
  *amn = i;
}

// impose periodic boundary conditions for node labels on ring lattice for Watts-Strogatz model
int wrap(int r, int n)
{
  if(r < 0) return n+r;
  if(r > n-1) return r-n;
  return r;
}

// build Watts-Strogatz-based networks matching n and e statistics
// there are several ways to do this.
void AMWattsStrogatz(int n, int e, int *amlist, int *amn, double beta)
{
  // the picture here is that we specify the mean degree k. we should choose this to give roughly e/n. this will usually not be an integer, so for every node we'll choose either the ceiling or the floor.
  double k = (double)e/n;
  int klo = (int)k;
  int khi = klo+1;
  double kd = k-klo;
  int thisk;
  int i, j, l;
  
  *amn = 0;
  
  for(i = 0; i < n; i++)
    {
      // decide whether to use the ceiling or the floor of mean degree for this node
      if(RND < kd) thisk = khi;
      else thisk = klo;

      // loop through this node's degree
      for(l = 1; l <= thisk; l++)
	{
	  // all edges start from this node
	  amlist[2*(*amn)+0] = i;
	  // edges going right may be rewired
	  if(RND < beta)
	    {
	      do{ j = RND*n; }while(j == i);
	      amlist[2*(*amn)+1] = j;
	    }
	  else
	    amlist[2*(*amn)+1] = wrap(i+l, n);
	  (*amn)++;
	}
    }
}


// structure and comparison function used in the quicksort codew for the geometric random graph
typedef struct { int a, b; double d; } Element;

int cmpfunc(const void *a, const void *b)
{
  Element A, B;
  A = *(Element*)a; B = *(Element*)b;
  return (A.d > B.d) - (A.d < B.d);
}

// build geometric random graph for given n, e
void AMGeometric(int n, int e, int *amlist, int *amn)
{
  double *x, *y;
  int i, j;
  Element *elements;
  int counter;
  double dx, dy;
  
  x = (double*)malloc(sizeof(double)*n);
  y = (double*)malloc(sizeof(double)*n);
  elements = (Element*)malloc(sizeof(Element)*n*n);

  (*amn) = 0;
  // place points randomly in unit square
  for(i = 0; i < n; i++)
    {
      x[i] = RND; y[i] = RND;
    }
  // compute distances between each pair, storing distance and pair labels in a structure
  counter = 0;
  for(i = 0; i < n; i++)
    {
      for(j = i+1; j < n; j++)
	{
	  // periodic boundary conditions
	  dx = (x[i]-x[j]); dy = (y[i]-y[j]);
	  if(dx > 0.5) dx = (x[i]-1)-x[j];
  	  if(dy > 0.5) dy = (y[i]-1)-y[j];
  	  if(dx < -0.5) dx = (x[i]+1)-x[j];
  	  if(dy < -0.5) dy = (y[i]+1)-y[j];

	  elements[counter].d = dx*dx+dy*dy;
	  elements[counter].a = i;
	  elements[counter].b = j;
	  counter++;
	}
    }
  // sort that structure by increasing distance
  qsort(elements, counter, sizeof(Element), cmpfunc);

  // go through the first e elements in the sorted structure, adding these pair labels to our adjacency matric
  for(i = 0; i < e; i++)
    {
      amlist[2*(*amn)+0] = elements[i].a;
      amlist[2*(*amn)+1] = elements[i].b;
      (*amn)++;
    }

  free(x);
  free(y);
  free(elements);
}
   
      

// build scale-free networks matching n and e statistics
// there are several ways to do this.
// flavour = 0 starts with n disconnected nodes and gives a node i a new edge with probability ~ deg(i)+1
// flavour = 1 starts with a linear graph and proceeds as with 0
// flavour = 2 does preferential attachment
void AMScaleFree(int n, int e, int *amlist, int *amn, int flavour)
{
  int *am;
  int i, j;
  int *deg;
  double *cumsum;
  double r;
  int nedge;
  
  // allocate memory for positions and convenient matrix
  deg = (int*)malloc(sizeof(int)*n);
  cumsum = (double*)malloc(sizeof(double)*n);
  am = (int*)malloc(sizeof(int)*n*n);
  
  // initialise empty encounter matrix and degrees
  for(i = 0; i < n; i++)
    deg[i] = 0;
  for(i = 0; i < n*n; i++)
    am[i] = 0;
  nedge = 0;

  // for flavour == 0, start with empty graph
  // for flavour == 1, start with a linear graph
  if(flavour == 1)
    {
      for(; nedge < n-1; nedge++)
	{
	  // add edge
	  amlist[2*nedge+0] = nedge;
	  amlist[2*nedge+1] = nedge+1;
	  am[nedge*n+(nedge+1)] = am[(nedge+1)*n+nedge] = 1;
	  deg[nedge]++; deg[nedge+1]++;
	}
    }
  // for flavour == 2, start by doing preferential attachment
  if(flavour == 2)
    {
      // assume we already have node 0
      // loop through the nodes to "add"
      for(i = 1; i < n; i++)
	{
	  // compute cumulative sum of degrees over "added" nodes
	  cumsum[0] = deg[0];
	  for(j = 1; j < i; j++)
	    cumsum[j] = cumsum[j-1]+deg[j]+1;
	  // roulette wheel selection for the node to connect our new node (i) to
	  r = RND*cumsum[i-1];
	  for(j = 0; r > cumsum[j]; j++);
	  // add edge
	  amlist[2*nedge+0] = i;
	  amlist[2*nedge+1] = j;
	  am[i*n+j] = am[j*n+i] = 1;
	  nedge++;
	  deg[i]++; deg[j]++;
	}
    }
  // keep adding edges up to e 
  for(; nedge < e; nedge++)
    {
      // compute cumulative sum of degrees over nodes
      cumsum[0] = deg[0];
      for(i = 1; i < n; i++)
	cumsum[i] = cumsum[i-1]+deg[i]+1;
      // roulette wheel selection for two nodes to connect
      r = RND*cumsum[n-1];
      for(i = 0; r > cumsum[i]; i++);
      do{
        r = RND*cumsum[n-1];
        for(j = 0; r > cumsum[j]; j++);
      }while(i == j || am[i*n+j] == 1);
      // add edge
      amlist[2*nedge+0] = i;
      amlist[2*nedge+1] = j;
      deg[i]++; deg[j]++;
      am[i*n+j] = am[j*n+i] = 0;
    }

  (*amn) = nedge;
  free(am);
  free(cumsum);
  free(deg);
}


// create an adjacency matrix matching given n, e statistics from a diffusive physical simulation
int AMFromSimulation(int n, int e, int *amlist, int *amn, double D, double kon, double koff, double V, double dmito, double kmito, int output, char *fname)
{
  // physical parameters -- somewhat arbitrary
  double cx = 100, cy = 30; // cell size
  double thresh = 1;        // encounter threshold

  double *x, *y;
  double *dx, *dy;
  int *active;
  double rhoon, rhooff;
  int *am;
  int i, j;
  long int timer;
  double *scale;
  FILE *fp;
  double r2;
  
  if(output != 0)
    {
      fp = fopen(fname, "w");
      fprintf(fp, "frame,mito,x,y\n");
    }

  if(INACTIVITY) { rhoon = 0.01; rhooff = 0.1; }
  else rhoon = rhooff = 0;
  
  timer = 0;
  // allocate memory for positions and convenient matrix
  x = (double*)malloc(sizeof(double)*n);
  y = (double*)malloc(sizeof(double)*n);
  dx = (double*)malloc(sizeof(double)*n);
  dy = (double*)malloc(sizeof(double)*n);
  active = (int*)malloc(sizeof(int)*n);
  am = (int*)malloc(sizeof(int)*n*n);
  scale = (double*)malloc(sizeof(double)*n);
  
  // initialise empty encounter matrix and random positions
  for(i = 0; i < n*n; i++)
    am[i] = 0;
  for(i = 0; i < n; i++)
    {
      x[i] = RND*cx;
      y[i] = RND*cy;
      dx[i] = dy[i] = 0;
      if(INACTIVITY)
        active[i] = (i < n/10);
      else
	active[i] = 1;
    }
  // loop while we have fewer than desired encounters
  *amn = 0;
  while(*amn < e && timer < 1e4)
    {
      timer++;

      for(i = 0; i < n; i++)
	scale[i] = 1;
      
      // look for encounters
      for(i = 0; i < n; i++)
	{
	  for(j = i+1; j < n; j++)
	    {
	      if(active[i] && active[j])
		{
		  r2 = (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]);
		  if(r2 < dmito*dmito)
		    {
		      scale[i] = kmito;
		      scale[j] = kmito;
		    }
		  if(*amn < e && am[i*n+j] == 0)
		    {
		      if(r2 < thresh*thresh)
			{
			  // populate encounter matrix and adj mat
			  am[i*n+j] = am[j*n+i] = 1;
			  amlist[2*(*amn)+0] = i;
			  amlist[2*(*amn)+1] = j;
			  (*amn)++;
			}
		    }
		}
	    }
	}

      //      if(timer % 1000 == 0)
      //	printf("%li\n", timer);
      // put a uniform diffusion kernel on each position
      for(i = 0; i < n; i++)
	{
	  if(active[i])
	    {
	      if(RND < rhooff) active[i] = 0;
	    }
	  else
	    {
	      if(RND < rhoon) active[i] = 1;
	    }
	  if(active[i])
	    {
	      // attach/detach from cytoskeleton
	      if(dx[i] == 0 && dy[i] == 0 && RND < kon)
		{
		  // choose horizontal or vertical motion randomly
		  if(RND < 0.5) { dx[i] = (RND < 0.5 ? V : -V); dy[i] = 0; }
		  else { dx[i] = 0; dy[i] = (RND < 0.5 ? V : -V); }
		}
	      else if((dx[i] != 0 || dy[i] != 0) && RND < koff)
		{
		  // back to diffusion
		  dx[i] = 0; dy[i] = 0;
		}

	      if(dx[i] == 0 && dy[i] == 0)
		{
		  x[i] += gsl_ran_gaussian(2.*D)*scale[i];
		  y[i] += gsl_ran_gaussian(2.*D)*scale[i];
		}
	      else
		{
		  x[i] += dx[i]*scale[i];
		  y[i] += dy[i]*scale[i];
		}
	      if(x[i] < 0) x[i] = 0;
	      if(x[i] > cx) x[i] = cx;
	      if(y[i] < 0) y[i] = 0;
	      if(y[i] > cy) y[i] = cy;
	    }
	}
      if(output != 0)
	{
	  for(i = 0; i < n; i++)
	    fprintf(fp, "%li,%i,%i,%.3f,%.3f\n", timer, i, active[i], x[i], y[i]);
	}
    }
  // free up memory
  free(x); free(y); free(dx); free(dy); free(am); free(scale); free(active);
  if(output != 0)
    fclose(fp);						    

  if(*amn < e) return -1;
  else return 0;
}						  


void mydebug(int i) { printf("wtf %i\n", i); exit(0); }

// play the bingo game itself, given a graph described by amlist with amn entries, max node label max, and number of genes L
// total and bingo store time series of total bingo score and number of bingos
void PlayBingo(int *amlist, int amn, int max, int L, double master, int bingotype, int *total, int *bingo)
{
  int i, j, k;
  int r;
  int *G, *H;
  int besttime;
  int t;
  int tmp;
  int score;
  int *edges;
  int nedges;
  
  edges = (int*)malloc(sizeof(int)*(amn+1));
  G = (int*)malloc(sizeof(int)*max*(L+1));
  H = (int*)malloc(sizeof(int)*max*(L+1));

  // populate initial node barcodes
  for(i = 0; i < max; i++)
    {
      if(RND > master)
	{
	  r = RND*L;
	  for(k = 0; k < L; k++)
	    {
	      G[i*L+k] = (k == r);
	      H[i*L+k] = (k == r);
	    }
	}
      else
	{
	  for(k = 0; k < L; k++)
	    {
	      G[i*L+k] = 1;
	      H[i*L+k] = 1;
	    }
	}
    }

  besttime = -1;

  for(i = 0; i < amn; i++)
    edges[i] = i;
  nedges = amn;
  
  // loop through number of exchange events
  for(t = 0; t < amn; t++)
    {
      // pick a random edge, or the earliest edge, from the remaining list
      if(bingotype == 1)
        r = RND*nedges;
      else
	r = 0;
      i = amlist[2*edges[r]+0];
      j = amlist[2*edges[r]+1];
      
      // pop this edge from the list
      for(k = r; k < nedges-1; k++)
	edges[k] = edges[k+1];
      nedges--;
      
      // give each node's history the other's genome
      for(k = 0; k < L; k++)
	{
	  if(G[i*L+k]) H[j*L+k] = 1;
	  if(G[j*L+k]) H[i*L+k] = 1;
	}

      // swap node genomes
      for(k = 0; k < L; k++)
	{
	  tmp = G[i*L+k];
	  G[i*L+k] = G[j*L+k];
	  G[j*L+k] = tmp;
	}

      // summary statistics -- highest and mean barcode occupancy
      total[t] = bingo[t] = 0;
      for(i = 0; i < max; i++)
	{
	  score = 0;
	  for(k = 0; k < L; k++)
	    score += H[i*L+k];
	  if(score == L)
	    bingo[t]++;
	  //	  if(bingo[0] > 0)
	  //  mydebug(4);
	  total[t] += score;
	}
    }

  free(edges);
  free(G);
  free(H);
}

void igraphStats(int *amlist, int amn, int n, int L, double *efficiency, int *cc, double *modularity, int *singletons, int *smalls, double *ccmean)
{
#ifdef _USEIGRAPH_
  igraph_t g;
  int i, j;
  double meand, norm;
  int numc;

  // A, B, C are sections for different network statistics
  
  // construct igraph graph from amlist
  igraph_empty(&g, n, 0);
  for(i = 0; i < amn; i++)
    igraph_add_edge(&g, amlist[2*i+0], amlist[2*i+1]);

  // A. global efficiency, using old igraph functions
  // first compute matrix (res) of distances between nodes
  // then manually compute efficiency
  
  igraph_matrix_t res;

  igraph_matrix_init(&res, 0, 0);
  igraph_shortest_paths_dijkstra(&g, &res, igraph_vss_all(), igraph_vss_all(), NULL, IGRAPH_OUT);

  meand = norm = 0;
  for(i = 0; i < n; i++)
    {
      for(j = i+1; j < n; j++)
	{
	  meand += 1./MATRIX(res, i, j);
	  norm++;
	}
    }
  *efficiency = meand/norm;
  
  // A. cleanup
  igraph_matrix_destroy(&res);

  // B. now, component structure. numc will store number of components. csize stores their sizes, which we loop through to get mean and number of singletons
  
  igraph_vector_t csize;

  igraph_vector_init(&csize, 0);
  igraph_clusters(&g, NULL, &csize, &numc, (igraph_connectedness_t) 0);
  *singletons = 0; *ccmean = 0; *smalls = 0;
  for(i = 0; i < igraph_vector_size(&csize); i++)
    {
      if(VECTOR(csize)[i] == 1) (*singletons)++;
      if(VECTOR(csize)[i] < L) (*smalls) += VECTOR(csize)[i];
      (*ccmean) += VECTOR(csize)[i];
    }
  (*ccmean) /= igraph_vector_size(&csize);
  *cc = numc;

  // B. cleanup
  igraph_vector_destroy(&csize);

  // C. community structure using walktrap. the modularity score is that coming from the final merge iteration of the algorithm
  
  igraph_matrix_t merges;
  igraph_vector_t tmpmodularity;
  long int no_of_nodes;
    
  igraph_vector_init(&tmpmodularity, 0);
  igraph_matrix_init(&merges, 0, 0);

  igraph_community_walktrap(&g, 0 /* no weights */,
			    4 /* steps */,
			    &merges, &tmpmodularity,
			    /* membership=*/ 0);

  no_of_nodes = igraph_vcount(&g);
  //    printf("Merges:\n");
  for (i = 0; i < igraph_matrix_nrow(&merges); i++) {
    *modularity = VECTOR(tmpmodularity)[i];
    /*  printf("%2.1li + %2.li -> %2.li (modularity %4.2f)\n",
	(long int)MATRIX(merges, i, 0),
	(long int)MATRIX(merges, i, 1),
	no_of_nodes + i,
	VECTOR(modularity)[i]);*/
  }

  // C. cleanup
  igraph_matrix_destroy(&merges);
  igraph_vector_destroy(&tmpmodularity);

  // overall cleanup
  igraph_destroy(&g);
#else
  *efficiency = 0;
  *cc = 0;
  *modularity = 0;
  *singletons = 0;
  *smalls = 0;
  *ccmean = 0;
#endif
}

double sd(double ss, double ss2, int n)
{
  return sqrt((ss2 - (ss*ss)/n) / (n - 1));
}
    
int main(int argc, char *argv[])
{
  int *amlist;
  int *degree;
  int amn;
  int *G, *H;
  int i, j, k;
  int r, r1, r2;
  int tmp;
  int *total, *bingo;
  int t;
  int it, sit;
  FILE *fp, *fp1, *fp2;
  char fstr[100], outstr[100];
  int nnodes;
  int found;
  int expt, exptref;
  int besttime;
  int nedge;
  int Lvals[] = {2, 3, 4, 5, 6, 7, 8, 9, 10};
  int L, Lindex;
  int cliquesize, cliquetype;
  double master;
  double D, kon, koff, V;
  int success;
  double score;
  int cc, singletons, smalls, maxdegree, mindegree;
  double modularity, meandegree, sddegree, efficiency, ccmean;
  int range, *rangeset;
  int bingotype = BINGOTYPE;
  int output, output_graphs;
  char filelabel[1000];
  
  double sum_meandegree, sum_meandegree_2, sum_sddegree, sum_sddegree_2, sum_rangedegree, sum_rangedegree_2, sum_efficiency, sum_efficiency_2, sum_singletons, sum_singletons_2, sum_ccmean, sum_ccmean_2, sum_modularity, sum_modularity_2, sum_smalls, sum_smalls_2, sum_bingoend[5], sum_bingoend_2[5], sum_cc, sum_cc_2;

  srand48(121);
  
  if(argc < 2)
    {
      printf("Need a file to analyse!");
      exit(0);
    }
  for(i = 2; i < argc; i++)
    {
      if(!strcmp(argv[i], "--output-simulation-traces\0")) OUTPUT_SIMULATION_TRACES = 1;
      if(!strcmp(argv[i], "--remove-singletons\0")) REMOVE_SINGLETONS = 1;
      if(!strcmp(argv[i], "--inactivity\0")) INACTIVITY = 1;
      if(!strcmp(argv[i], "--output-trajectories\0"))
	{
	  if(i == argc-1) printf("Number of trajectories to output not specified!");
	  else
	    {
	      OUTPUT_TRAJECTORIES = atoi(argv[i+1]);
	      i++;
	    }
	}
    }

  printf("Attempting to use file %s\n", argv[1]);
  printf("Flags:\n  output-simulation-traces %i\n  remove-singletons %i\n  output-trajectories %i\n  inactivity %i\n\n", OUTPUT_SIMULATION_TRACES, REMOVE_SINGLETONS, OUTPUT_TRAJECTORIES, INACTIVITY);
  sprintf(filelabel, "%s-%i-%i-%i", argv[1], bingotype, REMOVE_SINGLETONS, INACTIVITY);
  
  // get minimum and maximum node labels
  //  StatsFromFile(argv[1], &min, &max);

  // allocate memory for adj mat and node barcodes
  amlist = (int*)malloc(sizeof(int)*2*MAXN*MAXN);
  degree = (int*)malloc(sizeof(int)*MAXN);
  rangeset = (int*)malloc(sizeof(int)*MAXN);

  total = (int*)malloc(sizeof(int)*MAXE);
  bingo = (int*)malloc(sizeof(int)*MAXE);

  // read in adj mat
  AMFromFile(argv[1], amlist, &amn, &nnodes);
  /*  for(i = 0; i < 2*amn; i++)
      amlist[i] -= min;
      max = max-min+1;*/

  nedge = amn;

  if(nnodes > MAXN || nedge > MAXE)
    {
      printf("Preallocated memory too small. Change MAXN and/or MAXE in the code.\n");
      exit(0);
    }

  // output stats
  printf("Found %i nodes and %i edges, ", nnodes, amn);
  singletons = 0;
  for(i = 0; i < nnodes; i++)
    {
      found = 0;
      for(j = 0; j < 2*amn; j++)
	{
	  if(amlist[j] == i) found = 1;
	}
      if(found == 0) singletons++;
    }
  printf("%i singletons, ", singletons);
  meandegree = 0;
  for(i = 0; i < nnodes; i++) degree[i] = 0;
  for(j = 0; j < 2*amn; j++) degree[amlist[j]]++;
  for(i = 0; i < nnodes; i++) meandegree += degree[i];
  printf("mean degree %.3f\n", meandegree/nnodes);

  if(OUTPUT_TRAJECTORIES > 0)
    {
      sprintf(fstr, "Output/%s-results-frame.txt", filelabel);
      fp = fopen(fstr, "w");
      fprintf(fp, "L,m,expt,it,n.edges,mean.score,num.bingos,prop.edges,prop.bingos\n");
    }
  
  sprintf(fstr, "Output/%s-results-overall.txt", filelabel);
  fp2 = fopen(fstr, "w");
  fprintf(fp2, "L,m,expt,mean.degree.mean,mean.degree.sd,sd.degree.mean,sd.degree.sd,range.degree.mean,range.degree.sd,efficiency.mean,efficiency.sd,modularity.mean,modularity.sd,singleton.count.mean,singleton.count.sd,small.count.mean,small.count.sd,num.cc.mean,num.cc.sd,mean.cc.size.mean,mean.cc.size.sd");
  for(i = 1; i <= 5; i++)
    fprintf(fp2, ",bingo.%i.mean,bingo.%i.sd", i, i);
  fprintf(fp2, "\n");
      
  // experiment is laid out as follows
  // loop over different L values (number of genetic fragments)
  // for each, loop through different proportions of master circles m
  // then loop through experiments (bio / different model graphs)
  for(Lindex = 0; Lindex <= 8; Lindex++) 
    {
      L = Lvals[Lindex];
      for(master = 0; master <= 0.02; master += 0.01)
	{
	  output_graphs = 0;
	  if(Lindex == 0 && master == 0)
	    {
	      sprintf(outstr, "Output/%s-graphs.txt", filelabel);
	      fp1 = fopen(outstr, "w");
	      output_graphs = 1;
	    }

	  for(expt = 0; expt <= 37; expt++)
	    {
	      sum_meandegree = sum_meandegree_2 = sum_sddegree = sum_sddegree_2 = sum_rangedegree = sum_rangedegree_2 = sum_efficiency = sum_efficiency_2 = sum_cc = sum_cc_2 = sum_modularity = sum_modularity_2 = sum_singletons = sum_singletons_2 = sum_smalls = sum_smalls_2 = sum_ccmean = sum_ccmean_2 = 0;
	      for(i = 0; i < 5; i++)
		sum_bingoend[i] = sum_bingoend_2[i] = 0;
		  
	      // loop through random instances
	      for(it = 0; it < NIT; it++)
		{
		  printf("%i %.3f %i %i\n", L, master, expt, it);

		  output = 0;
		  if(OUTPUT_SIMULATION_TRACES)
		    {
		      if(Lindex == 0 && master == 0 && it == 0)
			{
			  output = 1;
			  sprintf(fstr, "Output/%s-sim-%i.csv", filelabel, expt);
			}
		    }
				  
		  // these either take the bio network (expt==0) or construct networks according to different protocols
		  switch(expt)
		    {
		      //		    case 0: StatsFromFile(argv[1], &min, &max); AMFromFile(argv[1], amlist, &amn, &max); for(i = 0; i < 2*amn; i++) amlist[i] -= min; max = max-min+1; break;
		    case 0: AMFromFile(argv[1], amlist, &amn, &nnodes); break;
		    case 1: AMFromSimulation(nnodes, nedge, amlist, &amn, 0.02, 0, 0, 0, 0, 0, output, fstr); break;        // just low diffusion
		    case 2: AMFromSimulation(nnodes, nedge, amlist, &amn, 0.1, 0, 0, 0, 0, 0, output, fstr); break;          // just mid diffusion 
		    case 3: AMFromSimulation(nnodes, nedge, amlist, &amn, 1, 0, 0, 0, 0, 0, output, fstr); break;         // just high diffusion
		    case 4: AMFromSimulation(nnodes, nedge, amlist, &amn, 0.02, 0.1, 0.1, 0.1, 0, 0, output, fstr); break;  // low diffusion, low cytoskeleton speed
		    case 5: AMFromSimulation(nnodes, nedge, amlist, &amn, 0.1, 0.1, 0.1, 1, 0, 0, output, fstr); break;      // mid diffusion, mid cytoskeleton speed
		    case 6: AMFromSimulation(nnodes, nedge, amlist, &amn, 1, 0.1, 0.1, 10, 0, 0, output, fstr); break;   // high diffusion, high cytoskeleton speed
		    case 7: AMFromSimulation(nnodes, nedge, amlist, &amn, 0.02, 0.5, 0.1, 0.1, 0, 0, output, fstr); break;   // low diffusion, low cytoskeleton speed, cytoskeleton favoured
		    case 8: AMFromSimulation(nnodes, nedge, amlist, &amn, 0.1, 0.5, 0.1, 1, 0, 0, output, fstr); break;    // mid diffusion, mid cytoskeleton speed, cytoskeleton favoured
		    case 9: AMFromSimulation(nnodes, nedge, amlist, &amn, 1, 0.5, 0.1, 10, 0, 0, output, fstr); break;   // high diffusion, high cytoskeleton speed, cytoskeleton favoured
      		    case 10: AMFromSimulation(nnodes, nedge, amlist, &amn, 0.02, 0.1, 0.1, 0.1, 3, 0.3, output, fstr); break;   // low diffusion, low cytoskeleton speed, cytoskeleton favoured, slowing
		    case 11: AMFromSimulation(nnodes, nedge, amlist, &amn, 0.1, 0.1, 0.1, 1, 3, 0.3, output, fstr); break;    // mid diffusion, mid cytoskeleton speed, cytoskeleton favoured
		    case 12: AMFromSimulation(nnodes, nedge, amlist, &amn, 1, 0.1, 0.1, 10, 3, 0.3, output, fstr); break;   // high diffusion, high cytoskeleton speed, cytoskeleton favoured
		    case 13: AMErdosRenyi(nnodes, nedge, amlist, &amn); break;                    // ER
		    case 14: AMScaleFree(nnodes, nedge, amlist, &amn, 0); break;                  // SF flavour 0 (start disconnected, new edge with probability ~ deg(i)+1)
		    case 15: AMScaleFree(nnodes, nedge, amlist, &amn, 1); break;                  // SF flavour 1 (start linear connected, then as 0)
		    case 16: AMScaleFree(nnodes, nedge, amlist, &amn, 2); break;                  // SF flavour 2 (preferential attachment)
		    case 17: AMWattsStrogatz(nnodes, nedge, amlist, &amn, 0.1); break;                  // WS beta = 0.1
		    case 18: AMWattsStrogatz(nnodes, nedge, amlist, &amn, 0.2); break;                  // WS beta = 0.2
		    case 19: AMWattsStrogatz(nnodes, nedge, amlist, &amn, 0.5); break;                  // WS beta = 0.5
		    case 20: AMGeometric(nnodes, nedge, amlist, &amn); break;                          // geometric random graph
		    case 21: AMStar(nnodes, nedge, amlist, &amn); break;                          // star graph
		    default: exptref = expt-22; cliquetype = 2+(exptref%2); cliquesize = 3+(exptref/2)*5; AMCliquey(nnodes, nedge, amlist, &amn, cliquesize, cliquetype); // cliquey graphs
		    }
		  // for reference, the "default" option above currently labels by clique size:  22,23 = 3; 24,25 = 8; 26,27 = 13; 28,29 = 18; 30,31 = 23; 32,33 = 28; 34,35 = 33; 36,37 = 38
		  
		  // output example graph structure
		  if(output_graphs == 1)
		    {
		      fprintf(fp1, "%i %i %i %i\n", expt, it, nnodes, nnodes);
		      for(j = 0; j < amn; j++)
			fprintf(fp1, "%i %i %i %i\n", expt, it, amlist[2*j+0]+1, amlist[2*j+1]+1);
		    }

		  // play bingo and output time series, and compute network statistics
		  PlayBingo(amlist, amn, nnodes, L, master, (expt <= 12 ? bingotype : 1), total, bingo);
		  igraphStats(amlist, amn, nnodes, L, &efficiency, &cc, &modularity, &singletons, &smalls, &ccmean);
		  meandegree = sddegree = 0;
		  mindegree = 9999; maxdegree = 0;
		  for(i = 0; i < nnodes; i++) degree[i] = rangeset[i] = 0;
		  for(j = 0; j < 2*amn; j++) degree[amlist[j]]++;
		  for(i = 0; i < nnodes; i++)
		    {
		      meandegree += degree[i];
		      rangeset[degree[i]] = 1;
		      if(degree[i] < mindegree) mindegree = degree[i];
		      if(degree[i] > maxdegree) maxdegree = degree[i];
		    }
		  meandegree /= nnodes;
		  for(i = 0; i < nnodes; i++)
		    sddegree += (degree[i]-meandegree)*(degree[i]-meandegree);
		  sddegree = sqrt(sddegree/(nnodes-1));
		  range = 0;
		  for(i = 0; i < nnodes; i++)
		    range += rangeset[i];

		  // output sample trajectory
		  if(it < OUTPUT_TRAJECTORIES)
		    {
		      for(t = 0; t < amn; t+=10)
			{
			  fprintf(fp, "%i,%.3f,%i,%i,%i,%f,%i,%f,%f\n", L, master, expt, it, t, (double)total[t]/nnodes, bingo[t], (double)t/amn, (double)bingo[t]/nnodes);
			}
		      fprintf(fp, "\n");
		    }

		  // update computational formula for summary statistics
		  sum_meandegree += meandegree; sum_meandegree_2 += meandegree*meandegree;
		  sum_sddegree += sddegree; sum_sddegree_2 += sddegree*sddegree;
		  sum_rangedegree += range; sum_rangedegree_2 += range*range;
		  sum_efficiency += efficiency;	sum_efficiency_2 += efficiency*efficiency;
		  sum_modularity += modularity; sum_modularity_2 += modularity*modularity;
		  sum_singletons += singletons; sum_singletons_2 += singletons*singletons;
		  sum_smalls += smalls; sum_smalls_2 += smalls*smalls;
		  sum_cc += cc; sum_cc_2 += cc*cc;
		  sum_ccmean += ccmean; sum_ccmean_2 += ccmean*ccmean;
		  sum_bingoend[0] += bingo[(int)(amn/100)]; sum_bingoend_2[0] += bingo[(int)(amn/100)]*bingo[(int)(amn/100)];
		  sum_bingoend[1] += bingo[(int)(amn/20)]; sum_bingoend_2[1] += bingo[(int)(amn/20)]*bingo[(int)(amn/20)];
		  sum_bingoend[2] += bingo[(int)(amn/10)]; sum_bingoend_2[2] += bingo[(int)(amn/10)]*bingo[(int)(amn/10)];
		  sum_bingoend[3] += bingo[(int)(amn/2)]; sum_bingoend_2[3] += bingo[(int)(amn/2)]*bingo[(int)(amn/2)];
		  sum_bingoend[4] += bingo[amn-1]; sum_bingoend_2[4] += bingo[amn-1]*bingo[amn-1];	  
		}
	      // output summary statistics
	      fprintf(fp2, "%i,%.3f,%i,", L, master, expt);
	      fprintf(fp2, "%f,%f,", sum_meandegree/NIT, sd(sum_meandegree, sum_meandegree_2, NIT));
	      fprintf(fp2, "%f,%f,", sum_sddegree/NIT, sd(sum_sddegree, sum_sddegree_2, NIT));
	      fprintf(fp2, "%f,%f,", sum_rangedegree/NIT, sd(sum_rangedegree, sum_rangedegree_2, NIT));
	      fprintf(fp2, "%f,%f,", sum_efficiency/NIT, sd(sum_efficiency, sum_efficiency_2, NIT));
	      fprintf(fp2, "%f,%f,", sum_modularity/NIT, sd(sum_modularity, sum_modularity_2, NIT));
	      fprintf(fp2, "%f,%f,", sum_singletons/NIT, sd(sum_singletons, sum_singletons_2, NIT));
	      fprintf(fp2, "%f,%f,", sum_smalls/NIT, sd(sum_smalls, sum_smalls_2, NIT));
	      fprintf(fp2, "%f,%f,", sum_cc/NIT, sd(sum_cc, sum_cc_2, NIT));
	      fprintf(fp2, "%f,%f", sum_ccmean/NIT, sd(sum_ccmean, sum_ccmean_2, NIT));
	      for(i = 0; i < 5; i++)
		fprintf(fp2, ",%f,%f", sum_bingoend[i]/NIT, sd(sum_bingoend[i], sum_bingoend_2[i], NIT));
	      fprintf(fp2, "\n");
	    }
	      

	  if(Lindex == 0 && master == 0)
	    {
	      fclose(fp1);
	    }
	}
    }
  fclose(fp);
  fclose(fp2);
      
  return 0;
}
       
