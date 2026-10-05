#ifndef GEOMETRY_C
#define GEOMETRY_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>

#include"../include/geometry.h"
#include"../include/gparam.h"

// the index functions declared in geometry.h (set by init_indexing_lexeo)
long (*cart_to_si)(int const * const cartcoord, GParam const * const param);
void (*si_to_cart)(int *cartcoord, long si, GParam const * const param);
long (*lex_to_si)(long lex, GParam const * const param);
long (*si_to_lex)(long si, GParam const * const param);
long (*sisp_and_t_to_si_compute)(long sisp, int t, GParam const * const param);
void (*si_to_sisp_and_t_compute)(long *sisp, int *t, long si, GParam const * const param);
long (*cart_to_si_rect)(int const * const cartcoord, Rectangle const * const most_update);

// single index 4d = even/odd lexicographic index 4d
// single index 3d = even/odd lexicographic index 3d
void init_indexing_lexeo(void)
  {
  cart_to_si = &cart_to_lexeo;
  si_to_cart = &lexeo_to_cart;
  lex_to_si = &lex_to_lexeo;
  si_to_lex = &lexeo_to_lex;
  sisp_and_t_to_si_compute=&lexeosp_and_t_to_lexeo;
  si_to_sisp_and_t_compute=&lexeo_to_lexeosp_and_t;
	
	// for rectangles
	cart_to_si_rect = &cart_to_lexeo_rect;
  }


// initialize geometry
void init_geometry(Geometry *geo, GParam const * const param)
  {
  int i, value, valuep, valuem, err;
  long r, rm, rp;
  int cartcoord[STDIM];

  // allocate memory
  err=posix_memalign((void**)&(geo->d_nnp), (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(long *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  err=posix_memalign((void**)&(geo->d_nnm), (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(long *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(param->d_volume); r++)
     {
     err=posix_memalign((void**)&(geo->d_nnp[r]), (size_t)INT_ALIGN, (size_t) STDIM * sizeof(long));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     err=posix_memalign((void**)&(geo->d_nnm[r]), (size_t)INT_ALIGN, (size_t) STDIM * sizeof(long));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  err=posix_memalign((void**)&(geo->d_timeslice), (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(int));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  err=posix_memalign((void**)&(geo->d_spacecomp), (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(long));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  err=posix_memalign((void**)&(geo->d_tsp), (size_t)INT_ALIGN, (size_t) param->d_size[0] * sizeof(long *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<param->d_size[0]; r++)
     {
     err=posix_memalign((void**)&(geo->d_tsp[r]), (size_t)INT_ALIGN, (size_t) param->d_space_vol * sizeof(long));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // INITIALIZE
  for(r=0; r<param->d_volume; r++)
     {
     si_to_cart(cartcoord, r, param);

     for(i=0; i<STDIM; i++)
        {
        value=cartcoord[i];

        valuep=value+1;
        if(valuep >= param->d_size[i])
          {
          valuep-=param->d_size[i];
          }
        cartcoord[i]=valuep;
        rp=cart_to_si(cartcoord, param);
        geo->d_nnp[r][i]=rp;

        valuem=value-1;
        if(valuem<0)
          {
          valuem+=param->d_size[i];
          }
        cartcoord[i]=valuem;
        rm=cart_to_si(cartcoord, param);
        geo->d_nnm[r][i]=rm;

        cartcoord[i]=value;
        }
     } // end of loop on r

  for(r=0; r<param->d_volume; r++)
     {
     si_to_sisp_and_t_compute(&rp, &value, r, param);
     geo->d_spacecomp[r]=rp;
     geo->d_timeslice[r]=value;
     geo->d_tsp[value][rp]=r;
     }

  #ifdef DEBUG
    test_geometry(geo, param);
  #endif
  }  


// free memory
void free_geometry(Geometry *geo, GParam const * const param)
  {
  long r;

  for(r=0; r<param->d_volume; r++)
     {
     free(geo->d_nnp[r]);
     free(geo->d_nnm[r]);
     }
  free(geo->d_nnp);
  free(geo->d_nnm);

  free(geo->d_timeslice);
  free(geo->d_spacecomp);
  for(r=0; r<param->d_size[0]; r++)
     {
     free(geo->d_tsp[r]);
     }
  free(geo->d_tsp);
  }


long nnp(Geometry const * const geo, long r, int i);


long nnm(Geometry const * const geo, long r, int i);


long sisp_and_t_to_si(Geometry const * const geo, long sisp, int t);


void si_to_sisp_and_t(long *sisp, int *t, Geometry const * const geo, long si);


void test_geometry(Geometry const * const geo, GParam const * const param)
  {
  long si, ris_test, si_bis, sisp;
  int dir, cart[STDIM], cartsp[STDIM-1], t;

  // test of lex_to_cart <-> cart_to_lex
  for(si=0; si < param->d_volume; si++)
     {
     lex_to_cart(cart, si, param);
     ris_test=cart_to_lex(cart, param);

     if(si != ris_test)
       {
       fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // test of lexeo_to_cart <-> cart_to_lexeo
  for(si=0; si < param->d_volume; si++)
     {
     lexeo_to_cart(cart, si, param);
     ris_test=cart_to_lexeo(cart, param);

     if(si != ris_test)
       {
       fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // test of nnp <-> nnm
  for(si=0; si < param->d_volume; si++)
     {
     for(dir=0; dir<STDIM; dir++)
        {
        si_bis=nnp(geo, si, dir);
        ris_test=nnm(geo, si_bis, dir);

        if(si != ris_test)
          {
          fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
        }
     }

  // test of lexsp_to_cartsp <-> cartsp_to_lexsp
  for(sisp=0; sisp < param->d_space_vol; sisp++)
     {
     lexsp_to_cartsp(cartsp, sisp, param);
     ris_test=cartsp_to_lexsp(cartsp, param);

     if(sisp != ris_test)
       {
       fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // test of lexeosp_to_cartsp <-> cartsp_to_lexeosp
  for(sisp=0; sisp < param->d_space_vol; sisp++)
     {
     lexeosp_to_cartsp(cartsp, sisp, param);
     ris_test=cartsp_to_lexeosp(cartsp, param);

     if(sisp != ris_test)
       {
       fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

   // test of lexeosp_and_t_to_lexeo <-> lexeo_to_lexeosp_and_t
   for(si=0; si<param->d_volume; si++)
      {
      lexeo_to_lexeosp_and_t(&sisp, &t, si, param);
      ris_test=lexeosp_and_t_to_lexeo(sisp, t, param);

      if(si != ris_test)
        {
        fprintf(stderr, "Problems while testing geometry! (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      }
  }



//------------ these are not to be used outside geometry.c ----------------

// cartesian coordinates -> lexicographic index
long cart_to_lex(int const * const cartcoord, GParam const * const param)
  {
  int i;
  long ris, aux;

  ris=0;
  aux=1;
  for(i=0; i<STDIM; i++)
     {
     ris+=cartcoord[i]*aux;
     aux*=param->d_size[i];
     }

  // ris = cartcoord[0]
  //      +cartcoord[1]*size[0]
  //      +cartcoord[2]*size[0]*size[1]
  //      +...
  //      +cartcoord[STDIM-1]*size[0]*size[1]*...*size[STDIM-2]

  return ris;

  /*
  int i;
  long ris, aux;

  ris=0;
  aux=1;
  for(i=STDIM-1; i>=0; i--)
     {
     ris+=cartcoord[i]*aux;
     aux*=param->d_size[i];
     }

  // ris = cartcoord[STDIM-1] +
  //      +cartcoord[STDIM-2]*size[STDIM-1]+
  //      +cartcoord[STDIM-3]*size[STDIM-1]*size[STDIM-2]+
  //      +...
  //      +cartcoord[0]*size[STDIM-1]*size[STDIM-2]*...*size[1]

  return ris;
  */
  }


// lexicographic index -> cartesian coordinates
void lex_to_cart(int *cartcoord, long lex, GParam const * const param)
  {
  int i;
  long aux[STDIM];

  aux[0]=1;
  for(i=1; i<STDIM; i++)
     {
     aux[i]=aux[i-1]*param->d_size[i-1];
     }
  // aux[0]=1
  // aux[1]=size[0]
  // aux[2]=size[0]*size[1]
  // ...
  // aux[STDIM-1]=size[0]*size[1]*...*size[STDIM-2]

  for(i=STDIM-1; i>=0; i--)
     {
     cartcoord[i]=(int) (lex/aux[i]);
     lex-=aux[i]*cartcoord[i];
     }

  /*
  int i;
  long aux[STDIM];

  aux[STDIM-1]=1;
  for(i=STDIM-2; i>=0; i--)
     {
     aux[i]=aux[i+1]*param->d_size[i+1];
     }
  // aux[STDIM-1] = 1
  // aux[STDIM-2] = size[STDIM-1]
  // aux[STDIM-3] = size[STDIM-1]*size[STDIM-2]
  // ...
  // aux[0]       = size[STDIM-1]*size[STDIM-2]*...*size[1]

  for(i=0; i<STDIM; i++)
     {
     cartcoord[i]=(int) (lex/aux[i]);
     lex-=aux[i]*cartcoord[i];
     }
  */
  }


// cartesian coordinates -> lexicographic eo index
long cart_to_lexeo(int const * const cartcoord, GParam const * const param)
  {
  long lex;
  int i, eo;

  lex=cart_to_lex(cartcoord, param);

  eo=0;
  for(i=0; i<STDIM; i++)
     {
     eo+=cartcoord[i];
     }

  if(eo % 2==0)
    {
    return lex/2;
    }
  else
    {
    return (lex + param->d_volume)/2;
    }
  // even sites are written first
  }


// lexicographic eo index -> cartesian coordinates
void lexeo_to_cart(int *cartcoord, long lexeo, GParam const * const param)
  {
  long lex;
  int i, eo;

  if(param->d_volume % 2 == 0)
    {
    if(lexeo < param->d_volume/2)
      {
      lex=2*lexeo;
      }
    else
      {
      lex=2*(lexeo-param->d_volume/2);
      }
    lex_to_cart(cartcoord, lex, param);

    eo=0;
    for(i=0; i<STDIM; i++)
       {
       eo+=cartcoord[i];
       }
    eo = eo % 2;

    // this is to take care of the case d_volume is even but not
    // all the lattice extents are even
    if( (eo == 0 && lexeo >= param->d_volume/2) ||
        (eo == 1 && lexeo < param->d_volume/2) )
      {
      lex+=1;
      lex_to_cart(cartcoord, lex, param);
      }
    }
  else
    {
    if(lexeo <= param->d_volume/2)
      {
      lex=2*lexeo;
      }
    else
      {
      lex=2*(lexeo-param->d_volume/2)-1;
      }
    lex_to_cart(cartcoord, lex, param);
    }
  }


//  lexicographic index -> lexicographic eo index
long lex_to_lexeo(long lex, GParam const * const param)
  {
  int cartcoord[STDIM];

  lex_to_cart(cartcoord, lex, param);

  return cart_to_lexeo(cartcoord, param);
  }


//  lexicographic eo index -> lexicographic index
long lexeo_to_lex(long lexeo, GParam const * const param)
  {
  int cartcoord[STDIM];

  lexeo_to_cart(cartcoord, lexeo, param);

  return cart_to_lex(cartcoord, param);
  }


// spatial cartesian coordinates -> spatial lexicographic index
long cartsp_to_lexsp(int const * const ccsp, GParam const * const param)
  {
  // the index for the spatial cartesian coord. goes from 0 to STDIM-2 hence ccsp[STDIM-1]
  // cc   = t x1 x2 ... x_{STDIM-1}
  // ccsp =   x1 x2     x_{STDIM-1}
  int i;
  long ris, aux;

  ris=0;
  aux=1;
  for(i=0; i<STDIM-1; i++)
     {
     ris+=ccsp[i]*aux;
     aux*=param->d_size[i+1];
     }

  // ris = ccsp[0]
  //      +ccsp[1]*size[1]
  //      +ccsp[2]*size[1]*size[2]
  //      +...
  //      +ccsp[STDIM-2]*size[1]*size[2]*...*size[STDIM-2]

  return ris;

  /*
  int i;
  long ris, aux;

  ris=0;
  aux=1;
  for(i=STDIM-2; i>=0; i--)
     {
     ris+=ccsp[i]*aux;
     aux*=param->d_size[i+1];
     }

  // ris = ccsp[STDIM-2] +
  //      +ccsp[STDIM-3]*size[STDIM-1]+
  //      +ccsp[STDIM-4]*size[STDIM-1]*size[STDIM-2]+
  //      +ccsp
  //      +ccsp[0]*size[STDIM-1]*size[STDIM-2]*...*size[2]

  return ris;
  */
  }


// spatial lexicographic index -> spatial cartesian coordinates
void lexsp_to_cartsp(int *ccsp, long lexsp, GParam const * const param)
  {
  // the index for the spatial cartesian coord. goes from 0 to STDIM-2 hence ccsp[STDIM-1]
  // cc   = t x1 x2 ... x_{STDIM-1}
  // ccsp =   x1 x2     x_{STDIM-1}

  int i;
  long aux[STDIM-1];

  aux[0]=1;
  for(i=1; i<STDIM-1; i++)
     {
     aux[i]=aux[i-1]*param->d_size[i];
     }
  // aux[0]=1
  // aux[1]=size[1]
  // aux[2]=size[1]*size[2]
  // ...
  // aux[STDIM-2]=size[1]*size[2]*...*size[STDIM-2]

  for(i=STDIM-2; i>=0; i--)
     {
     ccsp[i]=(int) (lexsp/aux[i]);
     lexsp-=aux[i]*ccsp[i];
     }

  /*
  int i;
  long aux[STDIM-1];

  aux[STDIM-2]=1;
  for(i=STDIM-3; i>=0; i--)
     {
     aux[i]=aux[i+1]*param->d_size[i+2];
     }
  // aux[STDIM-2] = 1
  // aux[STDIM-3] = size[STDIM-1]
  // aux[STDIM-4] = size[STDIM-1]*size[STDIM-2]
  // ...
  // aux[0]       = size[STDIM-1]*size[STDIM-2]*...*size[2]

  for(i=0; i<STDIM-1; i++)
     {
     ccsp[i]=(int) (lexsp/aux[i]);
     lexsp-=aux[i]*ccsp[i];
     }
  */
  }


// spatial cartesian coordinates -> spatial lexicographic eo index
long cartsp_to_lexeosp(int const * const ccsp, GParam const * const param)
  {
  long lexsp;
  int i, eo;

  lexsp=cartsp_to_lexsp(ccsp, param);

  eo=0;
  for(i=0; i<STDIM-1; i++)
     {
     eo+=ccsp[i];
     }

  if(eo % 2==0)
    {
    return lexsp/2;
    }
  else
    {
    return (lexsp + param->d_space_vol)/2;
    }
  }


// spatial lexicographic eo index -> spatial cartesian coordinates
void lexeosp_to_cartsp(int *ccsp, long lexeosp, GParam const * const param)
  {
  long lexsp;
  int i, eo;

  if(param->d_space_vol % 2 == 0)
    {
    if(lexeosp < param->d_space_vol/2)
      {
      lexsp=2*lexeosp;
      }
    else
      {
      lexsp=2*(lexeosp - param->d_space_vol/2);
      }
    lexsp_to_cartsp(ccsp, lexsp, param);

    eo=0;
    for(i=0; i<STDIM-1; i++)
       {
       eo+=ccsp[i];
       }
    eo = eo % 2;

    if( (eo == 0 && lexeosp >= param->d_space_vol/2) ||
        (eo == 1 && lexeosp < param->d_space_vol/2) )
      {
      lexsp+=1;
      lexsp_to_cartsp(ccsp, lexsp, param);
      }
    }
  else
    {
    if(lexeosp <= param->d_space_vol/2)
      {
      lexsp=2*lexeosp;
      }
    else
      {
      lexsp=2*(lexeosp - param->d_space_vol/2)-1;
      }
    lexsp_to_cartsp(ccsp, lexsp, param);
    }
  }


// spatial lexicographic index -> spatial lexicographic eo index
long lexsp_to_lexeosp(long lexsp, GParam const * const param)
  {
  int ccsp[STDIM];

  lexsp_to_cartsp(ccsp, lexsp, param);

  return cartsp_to_lexeosp(ccsp, param);
  }


//  spatial lexicographic eo index -> spatial lexicographic index
long lexeosp_to_lexsp(long lexeosp, GParam const * const param)
  {
  int ccsp[STDIM];

  lexeosp_to_cartsp(ccsp, lexeosp, param);

  return cartsp_to_lexsp(ccsp, param);
  }


// lexicographic eo spatial and time -> lexicographic eo index
long lexeosp_and_t_to_lexeo(long lexeosp, int t, GParam const * const param)
  {
  int cc[STDIM];

  lexeosp_to_cartsp(cc+1, lexeosp, param);
  cc[0]=t;

  return cart_to_lexeo(cc, param);
  }


// lexicographic eo index -> lexicographic eo spatial and time
void lexeo_to_lexeosp_and_t(long *lexeosp, int *t, long lexeo, GParam const * const param)
  {
  int i, cc[STDIM], ccsp[STDIM-1];

  lexeo_to_cart(cc, lexeo, param);

  *t=cc[0];

  for(i=0; i<STDIM-1; i++)
     {
     ccsp[i]=cc[i+1];
     }

  *lexeosp=cartsp_to_lexeosp(ccsp, param);
  }
	
// 4d only geometry for rectangles to be used for hierarchical update during parallel tempering

// reduce a generic integer component in the interval [0,L_max-1]
int periodic_condition(int const coord, int const L_max)
  {
  if(coord<0) return (L_max+coord%L_max);
  else if(coord>=L_max) return (coord%L_max);
  else return coord;
  }

// cartesian -> lexicographic eo (on a given rectangle)
// this function should not be used outside geometry.c
long cart_to_lexeo_rect(int const * const cartcoord, Rectangle const * const most_update)
  {
		
  int i, eo=0;
  long ris=0, aux=1;

  for(i=0; i<STDIM; i++)
    {
    ris+=cartcoord[i]*aux;
    aux*=most_update->d_size_rect[i];
		eo+=cartcoord[i];
    }
	eo=eo%2;
	return (eo*(most_update->d_vol_rect)+ris)/2; // even sites first
  }

void init_rect(Rectangle *most_update, int const L_R, GParam const * const param)
  {
	// sizes of rectangle and ranges of rect coordinates
	int aux_L[STDIM], size_min[STDIM], size_max[STDIM];
	int i,err;
	
	// for each value of defect_dir, determine the three orthogonal directions to it
	int perp_dir[4][3] = { {1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2} };

	// rectangle sizes
	aux_L[param->d_defect_dir]=2*L_R+1;
	aux_L[perp_dir[param->d_defect_dir][0]]=2*L_R+(param->d_L_defect[0]);
	aux_L[perp_dir[param->d_defect_dir][1]]=2*L_R+(param->d_L_defect[1]);
	aux_L[perp_dir[param->d_defect_dir][2]]=2*L_R+(param->d_L_defect[2]);

	// min and max of rect coordinates along every direction
	size_min[param->d_defect_dir]=(param->d_size[param->d_defect_dir])-1-L_R;
	size_min[perp_dir[param->d_defect_dir][0]]=-L_R;
	size_min[perp_dir[param->d_defect_dir][1]]=-L_R;
	size_min[perp_dir[param->d_defect_dir][2]]=-L_R;

	size_max[param->d_defect_dir]=(param->d_size[param->d_defect_dir])+L_R;
	size_max[perp_dir[param->d_defect_dir][0]]=(param->d_L_defect[0])+L_R;
	size_max[perp_dir[param->d_defect_dir][1]]=(param->d_L_defect[1])+L_R;
	size_max[perp_dir[param->d_defect_dir][2]]=(param->d_L_defect[2])+L_R;

  // d_size_rect[i] must not exceed d_size[i]
  // if so, the exceeding dimension of the rectangle is just the respective dimension of the whole lattice
  // and the i-th coordinate just ranges from 0 to d_size[i]-1
	for(i=0;i<STDIM;i++)
		{
		if(aux_L[i]>=param->d_size[i])
			{
			aux_L[i]=param->d_size[i];
			size_min[i]=0;
			size_max[i]=param->d_size[i];
			}
		}
		
  // volume of rectangle
  long V=1;
	for(i=0;i<STDIM;i++)
		V*=aux_L[i];

  // allocate rectangle
  err=posix_memalign((void **) &(most_update->rect_sites), (size_t) INT_ALIGN, (size_t) V * sizeof(long));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the rectangle! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  // save dimensions and volume of the rectangle
  for(i=0;i<STDIM;i++)
		most_update->d_size_rect[i]=aux_L[i];
  most_update->d_vol_rect=V;
  
  // save lexeo index of sites of the rectangle
  int coord[STDIM];        // cartesian coordinates on the whole lattice
  int real_coord[STDIM];   // cartesian coordinates after using periodic conditions
  long r,r_rect;           // lexeo index on the rectangle and lexeo index on the whole lattice
  int coord_rect[STDIM];   // cartesian coordinates on the rectangle
	
  coord_rect[0]=0;
  for(coord[0]=size_min[0];coord[0]<size_max[0];coord[0]++)
		{
		coord_rect[1]=0;
		for(coord[1]=size_min[1];coord[1]<size_max[1];coord[1]++)
			{
			coord_rect[2]=0;
			for(coord[2]=size_min[2];coord[2]<size_max[2];coord[2]++)
				{
				coord_rect[3]=0;
				for(coord[3]=size_min[3];coord[3]<size_max[3];coord[3]++)
					{
					// compute the real coordinates on the periodic lattice
					for(i=0;i<STDIM;i++) 
						real_coord[i]=periodic_condition(coord[i],param->d_size[i]);
					
					r_rect=cart_to_si_rect(coord_rect,most_update); // from cartesian coordinates on the rectangle to lexeo index on the rectangle
					r=cart_to_si(real_coord,param);                 // from cartesian coordinates on the lattice to lexeo index on the lattice
					most_update->rect_sites[r_rect]=r;

					coord_rect[3]++;
		      } // end of z loop
				coord_rect[2]++;
				} // end of y loop
			coord_rect[1]++;
			} // end of x loop   
		coord_rect[0]++;
		} // end of t loop
  }

void free_rect(Rectangle *most_update)
  {
  free(most_update->rect_sites);
  }
	
void init_rect_hierarc(Rectangle **most_update, Rectangle **clover_rect, GParam const * const param)
	{
	if(param->d_N_hierarc_levels==0)
		{
		most_update=NULL;
		clover_rect=NULL;
		}
	else
		{
		int i,err;
		err=posix_memalign( (void **) most_update, (size_t) INT_ALIGN, (size_t) param->d_N_hierarc_levels * sizeof(Rectangle) );
		if(err!=0)
			{
			fprintf(stderr, "Problems in allocating the rectangle! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
			}
		for(i=0; i<param->d_N_hierarc_levels; i++)
			{
			init_rect(&((*most_update)[i]), param->d_L_rect[i], param);
			}
	
		#ifdef THETA_MODE
			err=posix_memalign( (void **) clover_rect, (size_t) INT_ALIGN, (size_t) param->d_N_hierarc_levels * sizeof(Rectangle) );
			if(err!=0)
				{
				fprintf(stderr, "Problems in allocating the rectangle! (%s, %d)\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
				}
			for(i=0; i<param->d_N_hierarc_levels; i++)
				init_rect(&((*clover_rect)[i]), 2+param->d_L_rect[i], param);
		#else
			clover_rect=NULL;
			(void) clover_rect; // just to suppress compiler warning of unused variable
		#endif
		}
	}

void free_rect_hierarc(Rectangle *most_update, Rectangle *clover_rect, GParam const * const param)
	{
	if(param->d_N_hierarc_levels==0)
		{
		// just to avoid compiler warning of unused variables
		(void) most_update;
		(void) clover_rect;
		}
	else
		{
		int i;
		for(i=0;i<param->d_N_hierarc_levels;i++)
			free_rect(&(most_update[i]));
		free(most_update);
		#ifdef THETA_MODE
			for(i=0;i<param->d_N_hierarc_levels;i++)
				free_rect(&(clover_rect[i]));
			free(clover_rect);
		#else
			(void) clover_rect; // just to avoid compiler warning of unused variable
		#endif
		}
	}

// compute the square distance between sites i and j on a periodic lattice (toroidal geometry)
double square_distance(long const i, long const j, GParam const * const param)
{
  int mu;
  int x[STDIM], y[STDIM];
  double d[STDIM], half_size[STDIM];
  double distance2=0.0;

  si_to_cart(x, i, param); // i --> x
  si_to_cart(y, j, param); // j --> y
  for (mu=0; mu<STDIM; mu++)
  {
    d[mu]=( (double) (abs(x[mu] - y[mu])) ); // d_mu = | x_mu - y_mu |
    half_size[mu] = ((double) param->d_size[mu])/2.0;
    if ( d[mu] > half_size[mu] ) d[mu] = ((double) param->d_size[mu]) - d[mu]; // periodic boundaries: if d_mu > L/2 ==> d_mu = L - d_mu
    distance2 += d[mu]*d[mu]; // distance^2 = sum_mu d_mu^2
  }
  return distance2;
}

#endif
