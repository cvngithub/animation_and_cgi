#include "StableFluidsSim.h"
#include <Eigen/LU>
void PRINT_S(char* X)
{
    std::cout << X << std::endl;
}
void PRINT_SCA(double& X)
{
    std::cout << X << std::endl;
}
void PRINT_VEC(VectorXs X)
{
    std::cout << X << std::endl << std::endl;
}
void PRINT_MAT2(MatrixXs X)
{
    std::cout << X << std::endl << std::endl;
}
void PRINT_IJ(int& i, int& j)
{
  std::cout <<  i  << "   " << j << "\n";
}
double LERP(const double& a,const double& b,const double& x)
{
  return (1.0 - x)*a+x*b; 
}
void StableFluidsSim::diffuseD(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt)
{
  assert((*x0 == *x0).all());
  
  scalar a = diff * dt * N * N;
  *x = *x0;
  //PRINT_S("StableFluidsSim::diffuseD start");
  for (int k = 0; k < 20; k++)
  {
    //scalar tmp;
    //scalar error = 0.0;
    for (int i = 1; i <= N; i++)
    {
      for (int j = 1; j <= N; j++) // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
      {
        // STUDENTS: You will certainly need code here, do diffuse for ((1, N), (1, N))
        //(*x)(i, j) = ((*x0)(i, j) + a * ((*x)(i, j+1) + (*x)(i, j - 1) + (*x)(i + 1, j) + (*x)(i - 1, j))) / (1 + 4 * a);
        (*x)(i, j) = ((*x0)(i, j) + a * ((*x0)(i-1, j) + (*x0)(i+1, j) + (*x0)(i, j-1) + (*x0)(i, j+1)))/(1+4*a);
        /*
        if((i == 1) && (j == 1))
          (*x)(i, j) += a * (                (*x0)(i+1, j)                 + (*x0)(i, j+1) - 2 * (*x0)(i,j));
        else if((i == 1) && (j == N))
          (*x)(i, j) += a * (                (*x0)(i+1, j) + (*x0)(i, j-1)                 - 2 * (*x0)(i,j));
        else if((i == N) && (j == 1)) 
          (*x)(i, j) += a * ((*x0)(i-1, j)                                 + (*x0)(i, j+1) - 2 * (*x0)(i,j));
        else if((i == N) && (j == N))
          (*x)(i, j) += a * ((*x0)(i-1, j)                 + (*x0)(i, j-1)                 - 2 * (*x0)(i,j));
        else if (i == 1)
          (*x)(i, j) += a * (                (*x0)(i+1, j) + (*x0)(i, j-1) + (*x0)(i, j+1) - 3 * (*x0)(i,j));
        else if (i == N)
          (*x)(i, j) += a * ((*x0)(i-1, j)                 + (*x0)(i, j-1) + (*x0)(i, j+1) - 3 * (*x0)(i,j));
        else if (j == 1)
          (*x)(i, j) += a * ((*x0)(i-1, j) + (*x0)(i+1, j)                 + (*x0)(i, j+1) - 3 * (*x0)(i,j));
        else if (j == N)
          (*x)(i, j) += a * ((*x0)(i-1, j) + (*x0)(i+1, j) + (*x0)(i, j-1)                 - 3 * (*x0)(i,j));
        else
          (*x)(i, j) += a * ((*x0)(i-1, j) + (*x0)(i+1, j) + (*x0)(i, j-1) + (*x0)(i, j+1) - 4 * (*x0)(i,j));
        */
      }
    }
    //PRINT_SCA(error);
  }
  //PRINT_S("StableFluidsSim::diffuseD end");
}

void StableFluidsSim::diffuseU(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt)
{
  assert((*x0 == *x0).all());
  //PRINT_S("StableFluidsSim::diffuseU start");
  scalar a = diff * dt * N * N;
  *x = *x0;
  for (int k = 0; k < 20; k++)
  {
    for (int i = 1; i <= N; i++)
    {
      for (int j = 0; j <= N; j++) // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
      {
        // STUDENTS: You will certainly need code here, do diffuse for ((1, N), (0, N)), note the case when (j == 0) or (j == N) need special treatment
          // q^n+1_i = a*(u((i,j-1) + 2u(i,j) + u(i,j+1)
        if (j == 0)
					(*x)(i, j) = ((*x0)(i, j) + a*((*x0)(i, j+1) + (*x0)(i+1, j) + (*x0)(i-1, j))) / (1+3*a);
				else if (j == N)
					(*x)(i, j) = ((*x0)(i, j) + a*((*x0)(i, j-1) + (*x0)(i+1, j) + (*x0)(i-1, j))) / (1+3*a);
				else
					(*x)(i, j) = ((*x0)(i, j) + a*((*x0)(i, j+1) + (*x0)(i, j-1) + (*x0)(i+1, j) + (*x)(i-1, j))) / (1+4*a);      
      }
    }
  }
  //PRINT_S("StableFluidsSim::diffuseU end");
}

void StableFluidsSim::diffuseV(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt)
{
  assert((*x0 == *x0).all());
  //PRINT_S("StableFluidsSim::diffuseV diffuseV");
  scalar a = diff * dt * N * N;
  *x = *x0;
  for (int k = 0; k < 20; k++)
  {
    for (int i = 0; i <= N; i++)
    {
      for (int j = 1; j <= N; j++) // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
      {       
        if(i == 0)
       		(*x)(i, j) = ((*x0)(i, j) + a*((*x0)(i, j+1) + (*x0)(i, j-1) + (*x0)(i+1, j))) / (1+3*a);
				else if (i == N)
					(*x)(i, j) = ((*x0)(i, j) + a*((*x0)(i, j+1) + (*x0)(i, j-1) + (*x0)(i-1, j))) / (1+3*a);
				else
					(*x)(i, j) = ((*x0)(i, j) + a*((*x0)(i, j+1) + (*x0)(i, j-1) + (*x0)(i+1, j) + (*x0)(i-1, j))) / (1+4*a);
        
      }
    }
  }
  //PRINT_S("StableFluidsSim::diffuseV end");
}

scalar StableFluidsSim::interpolateD(ArrayXs * d, scalar i, scalar j)
{
    // STUDENTS: You will certainly need code here, note the indices should be CLAMP-ed to (0, m_N), since we have to use (i + 1) and (j + 1)
  //PRINT_S("StableFluidsSim::diffuseD start");
  int i1 = CLAMP((int)(i-0.5), 0, m_N);//choose i1 as btw (1,N)
  int i2 = i1 + 1;
  int j1 = CLAMP((int)(j-0.5), 0, m_N);//choose j1 as btw (1,N)
  int j2 = j1 + 1;
  double s = CLAMP(i - i1, 0, 1);
  double t = CLAMP(j - j1, 0, 1);
  double d1 = LERP((*d)(i1 ,j1) , (*d)(i2 ,j1) , s);
  double d2 = LERP((*d)(i1 ,j2) , (*d)(i2 ,j2) , s);
  //PRINT_S("StableFluidsSim::diffuseD end");
  return LERP(d1 , d2 , t);
}

scalar StableFluidsSim::interpolateU(ArrayXs * u, scalar i, scalar j)
{
    // STUDENTS: You will certainly need code here,
    // note the i index should be CLAMP-ed to (0, m_N), 
    // while j index should be CLAMP-ed to (0, m_N-1), 
    // since we have to use (i + 1) and (j + 1)
  /* Document p.5
    3: i1 ← CLAMP(integer(i)) ∈ (0, N )
    4: i2←i1+1
    5: j1 ← CLAMP(integer(j − 0.5)) ∈ (0, N − 1)
    6: j2←j2+1
    7: s←CLAMP(i−i1)∈(0,1)
    8: t←CLAMP(j−j1 −0.5)∈(0,1)
    9: u1 ← LERP(ui1 ,j1 , ui2 ,j1 , s)
    10: u2 ← LERP(ui1 ,j2 , ui2 ,j2 , s)
    11: return LERP(u1 , u2 , t)
  */
  //PRINT_S("StableFluidsSim::interpolateU start");
  int i1 = CLAMP((int)i, 0, m_N);//choose i1 as btw (1,N)
  int i2 = i1 + 1;
  int j1 = CLAMP((int)(j-0.5), 0, m_N - 1);//choose j1 as btw (1,N-1)
  int j2 = j1 + 1;

  double s = CLAMP(i - i1, 0, 1);
  double t = CLAMP(j - j1 - 0.5, 0, 1);
  double u1 = LERP((*u)(i1 ,j1) , (*u)(i2 ,j1) , s);
  double u2 = LERP((*u)(i1 ,j2) , (*u)(i2 ,j2) , s);
  //PRINT_S("StableFluidsSim::interpolateU end");
  //PRINT_SCA(t);
  return LERP(u1 , u2 , t);
}
 
scalar StableFluidsSim::interpolateV(ArrayXs * v, scalar i, scalar j)
{
  // STUDENTS: You will certainly need code here
  //PRINT_S("StableFluidsSim::interpolateV start");
  //PRINT_S("-2");
  int i1 = CLAMP((int)(i-0.5), 0, m_N - 1);//choose i1 as btw 1,N-1)
  int i2 = i1 + 1;
  int j1 = CLAMP((int)(j), 0, m_N);//choose j1 as btw (1,N)
  int j2 = j1 + 1;

  double s = CLAMP(i - i1 - 0.5, 0, 1);
  double t = CLAMP(j - j1, 0, 1);
  double v1 = LERP((*v)(i1 ,j1) , (*v)(i1 ,j2) , t);
  double v2 = LERP((*v)(i2 ,j1) , (*v)(i2 ,j2) , t);
  //PRINT_S("StableFluidsSim::interpolateV end");
  //PRINT_SCA(t);
  return LERP(v1 , v2 , s);
}

void StableFluidsSim::advectU(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt)
{
  assert((*x0 == *x0).all());
  assert((*u == *u).all());
  assert((*v == *v).all());
  
  *x = *x0;
  //PRINT_S("StableFluidsSim::advectU start");
  for (int i = 1; i <= N; i++)
  {
    for (int j = 0; j <= N; j++)
    {
      // STUDENTS: You will certainly need code here,
      // formula for u stated on 8.3
      scalar dx = 0.5, dy = 0.0;
      //PRINT_S("StableFluidsSim::advectU");
      //PRINT_IJ(i, j);
      scalar bi =i + dy - dt * N * interpolateV(v, i + dy, j + dx);
      scalar bj =j + dx - dt * N * interpolateU(u, i + dy, j + dx);
      (*x)(i,j) = interpolateU(x0,bi,bj);
      
      // now you have the backward-traced velocity, minus it from the current position (i + 0, j + 0.5), then sample the velocity again.
    }
  }
  //PRINT_S("StableFluidsSim::advectU end");
}

void StableFluidsSim::advectV(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt)
{
    assert((*x0 == *x0).all());
    assert((*u == *u).all());
    assert((*v == *v).all());
  //PRINT_S("StableFluidsSim::advectV start");
  for (int i = 0; i <= N; i++)
  {
    for (int j = 1; j <= N; j++)
    {
      // STUDENTS: You will certainly need code here
      // stated on 8.3
      scalar dx = 0.0, dy = 0.5;
      //PRINT_S("StableFluidsSim::advectV");
      //PRINT_IJ(i, j);
      scalar bi =i + dy - dt * N * interpolateV(v, i + dy, j + dx);// whether swap?
      scalar bj =j + dx - dt * N * interpolateU(u, i + dy, j + dx);
      (*x)(i,j) = interpolateV(x0,bi,bj); 
    }
  }
  //PRINT_S("StableFluidsSim::advectV end");
}

void StableFluidsSim::advectD(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt)
{
  assert((*x0 == *x0).all());
  assert((*u == *u).all());
  assert((*v == *v).all());
  
  // STUDENTS: You will certainly need code here, advect for ((1, N), (1, N))
  //PRINT_S("StableFluidsSim::advectD start");
  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j <= N; j++)
    {
      // STUDENTS: You will certainly need code here
      // stated on 8.3
      scalar dx = 0.0, dy = 0.0;
      //PRINT_S("StableFluidsSim::advectD");
      //PRINT_IJ(i, j);
      scalar bi =i + dy - dt * N * interpolateV(v, i + dy, j + dx);
      scalar bj =j + dx - dt * N * interpolateU(u, i + dy, j + dx);
      (*x)(i,j) = interpolateD(x0,bi,bj);
    }
  }
  //PRINT_S("StableFluidsSim::advectD end");
}

void StableFluidsSim::project(int N, ArrayXs * u, ArrayXs * v, ArrayXs * u0, ArrayXs * v0)
{
  if (VERBOSE) std::cout << "u0: " << std::endl << *u0 << std::endl << std::endl;
  if (VERBOSE) std::cout << "v0: " << std::endl << *v0 << std::endl << std::endl;

  ArrayXs div(N + 2, N + 2);
  ArrayXs p(N + 2, N + 2);
  div.setZero();
  p.setZero();
  scalar h = 1.0 / N;
  
  // STUDENTS: You will certainly need code here
  // set solid boundary conditions, 
  // 0 the most top and bottom row 
  // left and right column of u0, v0
  //PRINT_S("StableFluidsSim::project BC");
  for (int i = 0; i <= N+1; i++)
  {
    // doc 5
    // ufluid = 0 for vertical walls
    // vfluid = 0 for horizontal walls
    (*u0)(i, 0) = 0.0;
    (*u0)(i, N) = 0.0;
    (*v0)(0, i) = 0.0;
    (*v0)(N, i) = 0.0;
  }
  for (int i = 0; i <= N; i++)
	{
		(*u0)(0, i) = 0;
		(*u0)(N+1, i) = 0;
		(*v0)(i, 0) = 0;
		(*v0)(i, N+1) = 0;
	}

  //PRINT_S("StableFluidsSim::project div");
  //int maxi, maxj;
  //double maxdiv = 0;
  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j <= N; j++)
    {
      // compute divergence of the velocity field, 
      // note the divergence field is available from ((1, N), (1, N))
      // formula 4.13 of fluides notes
      // div u(i,j) = u(i+0.5,j)-u(i-0.5,j))/dx + u(i,j+0.5)-u(i,j-0.5))/dx
      div(i, j) = ((*u0)(i,j) - (*u0)(i,j-1)) / h + ((*v0)(i,j) - (*v0)(i-1,j)) / h ;
      //if (maxdiv < div(i, j)){
      //  maxdiv = div(i, j);
      //  maxi = i;
      //  maxj = j;
      //}
    }
  }
  //PRINT_IJ(maxi, maxj);
  //PRINT_SCA(maxdiv);
  //PRINT_S("StableFluidsSim::project p");
  for (int k = 0; k < 20; k++)
  {
    //scalar tmp;
    //scalar error = 0.0;
    for (int i = 1; i <= N; i++)
    {
      for (int j = 1; j <= N; j++) // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
      {
        //tmp = p(i, j);
        // solve for pressure inside the region ((1, N), (1, N))
        // from the formula of doc section 5.
        // fig.5 pressure gradient across walls with yellow circles should also equal zero
        if ((i == 1) && (j == 1))
          p(i, j) = (p(i+1,j)            + p(i,j+1)            - h*h*div(i,j)) / 2.0;
        else if ((i == N) && (j == 1))
          p(i, j) = (           p(i-1,j) + p(i,j+1)            - h*h*div(i,j)) / 2.0;  
        else if ((i == 1) && (j == N))
          p(i, j) = (p(i+1,j)                       + p(i,j-1) - h*h*div(i,j)) / 2.0;
        else if ((i == N) && (j == N))
          p(i, j) = (           p(i-1,j)            + p(i,j-1) - h*h*div(i,j)) / 2.0;
        else if (j == 1)
          p(i, j) = (p(i+1,j) + p(i-1,j) + p(i,j+1)            - h*h*div(i,j)) / 3.0;
        else if (j == N)
          p(i, j) = (p(i+1,j) + p(i-1,j)            + p(i,j-1) - h*h*div(i,j)) / 3.0;
        else if (i == 1) 
          p(i, j) = (p(i+1,j)            + p(i,j+1) + p(i,j-1) - h*h*div(i,j)) / 3.0;
        else if (i == N)
          p(i, j) = (           p(i-1,j) + p(i,j+1) + p(i,j-1) - h*h*div(i,j)) / 3.0;
        else
          p(i, j) = (p(i+1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1) - h*h*div(i,j)) / 4.0;
        //error += (p(i, j) - tmp) * (p(i, j) - tmp);
        
      }
    }
    //PRINT_SCA(error);
  }
 

  (*u) = (*u0); 
  (*v) = (*v0);
  //PRINT_S("StableFluidsSim::project uv");
  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j < N; j++)
    {
      // apply pressure to correct velocities 
      // ((1, N), (1, N)) for u, ((1, N), (1, N)) for v
      // from the 2nd formula of doc section 5.
      (*u)(i, j)  = (*u0)(i, j) - (p(i, j+1) - p(i, j))/h;
                    //+ m_visc * ((*u0)(i, j+1) - 2.0 * (*u0)(i, j) + (*u0)(i, j-1)) / h/h;
  		(*v)(j, i) =  (*v0)(j, i) - (p(j+1, i) - p(j, i))/h;
                   // + m_visc * ((*v0)(j+1, i) - 2.0 * (*v0)(j, i) + (*v0)(j-1, i)) / h/h;
    }
  }
}

