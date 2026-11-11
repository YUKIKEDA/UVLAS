#include <cstdint>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif
using Clock=std::chrono::steady_clock;
double seconds(Clock::time_point a){return std::chrono::duration<double>(Clock::now()-a).count();}
std::uint64_t peak_bytes(){
#ifdef _WIN32
  PROCESS_MEMORY_COUNTERS p{};p.cb=sizeof(p);if(!GetProcessMemoryInfo(GetCurrentProcess(),&p,sizeof(p)))throw std::runtime_error("memory query");return p.PeakWorkingSetSize;
#else
  rusage p{};if(getrusage(RUSAGE_SELF,&p))throw std::runtime_error("memory query");return std::uint64_t(p.ru_maxrss)*1024;
#endif
}
int main(){try{
  Eigen::setNbThreads(1);std::cout.precision(17);
  using Mat=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>;
  constexpr Eigen::Index n=5000;auto start=Clock::now();Mat A(n,n);
  for(Eigen::Index j=0;j<n;++j)for(Eigen::Index i=0;i<n;++i)A(i,j)=(i==j?2.:0.)+.1/(1.+std::abs(i-j))+.00001*std::sin(.03*i+.07*j);
  Eigen::MatrixXd truth(n,8);for(Eigen::Index j=0;j<8;++j)for(Eigen::Index i=0;i<n;++i)truth(i,j)=std::sin(.01*i+j);
  Eigen::MatrixXd b=A*truth;double assembly=seconds(start);start=Clock::now();Eigen::PartialPivLU<Mat> lu(A);double factor=seconds(start);start=Clock::now();Eigen::MatrixXd x=lu.solve(b);double solve=seconds(start);
  double error=(x-truth).norm()/truth.norm(),residual=(A*x-b).norm()/b.norm();if(error>1e-11||residual>1e-11)throw std::runtime_error("dense capacity correctness");
  std::cout<<"{\"case\":\"dense_5000_8rhs_capacity_only\",\"assembly_seconds\":"<<assembly<<",\"factor_seconds\":"<<factor<<",\"solve_seconds\":"<<solve<<",\"relative_error\":"<<error<<",\"residual\":"<<residual<<"}"<<std::endl;
  // Genuine two-dimensional fill, still NOT the PG patch graph or a product performance case.
  constexpr std::int64_t nx=400,ny=250,nd=nx*ny;using Sp=Eigen::SparseMatrix<double,Eigen::ColMajor,std::int64_t>;
  start=Clock::now();std::vector<Eigen::Triplet<double,std::int64_t>> t;t.reserve(25*nd);
  for(std::int64_t y=0;y<ny;++y)for(std::int64_t z=0;z<nx;++z){auto i=y*nx+z;double diagonal=1.;for(int dy=-2;dy<=2;++dy)for(int dz=-2;dz<=2;++dz){if((!dy&&!dz)||y+dy<0||y+dy>=ny||z+dz<0||z+dz>=nx)continue;auto j=(y+dy)*nx+z+dz;double w=.2/(std::abs(dy)+std::abs(dz));t.emplace_back(i,j,-w);diagonal+=w;}t.emplace_back(i,i,diagonal);}
  Sp S(nd,nd);S.setFromTriplets(t.begin(),t.end());S.makeCompressed();std::vector<Eigen::Triplet<double,std::int64_t>>().swap(t);Eigen::VectorXd xt=Eigen::VectorXd::LinSpaced(nd,.1,1.),rhs=S*xt;
  double sparse_assembly=seconds(start);start=Clock::now();Eigen::SparseLU<Sp,Eigen::COLAMDOrdering<std::int64_t>> sl;sl.analyzePattern(S);sl.factorize(S);if(sl.info()!=Eigen::Success)throw std::runtime_error("sparse capacity factor");double sf=seconds(start);start=Clock::now();Eigen::VectorXd xs=sl.solve(rhs);double ss=seconds(start);double se=(xs-xt).norm()/xt.norm(),sr=(S*xs-rhs).norm()/rhs.norm();if(se>1e-11||sr>1e-11)throw std::runtime_error("sparse capacity correctness");
  std::cout<<"{\"case\":\"sparse_100k_2d_radius2_capacity_only\",\"nnz\":"<<S.nonZeros()<<",\"factor_nnz\":"<<sl.nnzL()+sl.nnzU()<<",\"assembly_seconds\":"<<sparse_assembly<<",\"factor_seconds\":"<<sf<<",\"solve_seconds\":"<<ss<<",\"relative_error\":"<<se<<",\"residual\":"<<sr<<",\"peak_process_bytes\":"<<peak_bytes()<<",\"passed\":true}"<<std::endl;return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
