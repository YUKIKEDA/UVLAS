#include <cstdint>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/DenseGenMatProd.h>
#include <expected>
#include <iostream>
#include <fstream>
#include <complex>
#include <limits>
#include <stdexcept>
#include <vector>
#include <cmath>
using Mat=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>;
using CM=Eigen::MatrixXcd;
using Sp=Eigen::SparseMatrix<double,Eigen::ColMajor,std::int64_t>;
void check(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
enum class Failure{Shape,Nonfinite,Singular};
std::expected<Mat,Failure> checked_solve(const Mat& A,const Mat& b){
  if(A.rows()!=A.cols()||A.rows()!=b.rows())return std::unexpected(Failure::Shape);
  if(!A.allFinite()||!b.allFinite())return std::unexpected(Failure::Nonfinite);
  if(A.rows()==0)return Mat(0,b.cols());
  Eigen::FullPivLU<Mat> lu(A);if(lu.rank()!=A.rows())return std::unexpected(Failure::Singular);
  Mat x=lu.solve(b);if(!x.allFinite())return std::unexpected(Failure::Nonfinite);return x;
}
struct Operator{
  using Scalar=double;const Mat& A;mutable std::int64_t calls=0;
  Eigen::Index rows()const{return A.rows();}Eigen::Index cols()const{return A.cols();}
  void perform_op(const double* x,double* y)const{++calls;Eigen::Map<Eigen::VectorXd> out(y,rows());out.noalias()=A*Eigen::Map<const Eigen::VectorXd>(x,cols());}
};
void test_arnoldi(const Mat& A,Eigen::Index nev,Eigen::Index ncv){
  Operator op{A};Spectra::GenEigsSolver<Operator> solver(op,nev,ncv);Eigen::VectorXd seed=Eigen::VectorXd::LinSpaced(A.rows(),.1,1.);solver.init(seed.data());auto count=solver.compute(Spectra::SortRule::LargestMagn,2000,1e-12);
  check(count==nev&&solver.info()==Spectra::CompInfo::Successful,"Arnoldi convergence");auto values=solver.eigenvalues();auto V=solver.eigenvectors();Eigen::EigenSolver<Mat> full(A,false);double maxres=0,maxdifference=0;
  for(Eigen::Index j=0;j<values.size();++j){double res=(A.cast<std::complex<double>>()*V.col(j)-values[j]*V.col(j)).norm()/(A.norm()*V.col(j).norm());maxres=std::max(maxres,res);double best=std::numeric_limits<double>::infinity();for(auto z:full.eigenvalues())best=std::min(best,std::abs(z-values[j]));maxdifference=std::max(maxdifference,best);}
  check(maxres<1e-10&&maxdifference<1e-8,"Arnoldi original residual/reference");std::cout<<"{\"case\":\"real_matrix_free_Arnoldi\",\"n\":"<<A.rows()<<",\"residual\":"<<maxres<<",\"full_eigen_difference\":"<<maxdifference<<",\"calls\":"<<op.calls<<"}\n";
}
int main(int argc,char** argv){try{
  static_assert(sizeof(Eigen::Index)==8);Eigen::setNbThreads(1);std::cout.precision(17);
  if(argc==2){std::ifstream in(argv[1]);Eigen::Index n=0;in>>n;check(in.good()&&n>=8&&n<=2048,"oracle size");Mat A(n,n);for(Eigen::Index i=0;i<n;++i)for(Eigen::Index j=0;j<n;++j)in>>A(i,j);check(in.good()&&A.allFinite(),"oracle data");test_arnoldi(A,6,std::min<Eigen::Index>(n,20));return 0;}
  check(checked_solve(Mat(0,0),Mat(0,2)).has_value(),"empty solve");check(!checked_solve(Mat::Zero(2,3),Mat::Zero(2,1)),"shape rejection");check(!checked_solve(Mat::Zero(2,2),Mat::Ones(2,1)),"singular rejection");Mat bad=Mat::Identity(2,2);bad(0,0)=std::numeric_limits<double>::quiet_NaN();check(!checked_solve(bad,Mat::Ones(2,1)),"NaN rejection");
  constexpr Eigen::Index n=48,r=4;Mat A(n,n);CM W(n,r),truth(n,3);for(Eigen::Index i=0;i<n;++i){for(Eigen::Index j=0;j<n;++j)A(i,j)=.02/(1.+std::abs(i-j))+.001*std::sin(.4*i+.7*j)+(i==j?2.:0.);for(Eigen::Index j=0;j<r;++j)W(i,j)={.03*std::sin(.2*i+j),.02*std::cos(i+.3*j)};for(Eigen::Index j=0;j<3;++j)truth(i,j)={std::sin(.3*(i+1)*(j+1)),std::cos(.2*i+.7*j)};}
  CM effective=A.cast<std::complex<double>>();effective.rightCols(r)+=W;CM rhs=effective*truth;Eigen::PartialPivLU<Mat> lu(A);CM Y=lu.solve(W.real()).cast<std::complex<double>>();Y+=std::complex<double>(0,1)*lu.solve(W.imag()).cast<std::complex<double>>();CM G=lu.solve(rhs.real()).cast<std::complex<double>>();G+=std::complex<double>(0,1)*lu.solve(rhs.imag()).cast<std::complex<double>>();CM small=CM::Identity(r,r)+Y.bottomRows(r);CM solved=G-Y*small.partialPivLu().solve(G.bottomRows(r));double schur_error=(solved-truth).norm()/truth.norm();check(schur_error<1e-12&&(effective*solved-rhs).norm()/rhs.norm()<1e-12,"complex rank update");std::cout<<"{\"case\":\"real_LU_complex_TE_update\",\"error\":"<<schur_error<<"}\n";
  constexpr Eigen::Index ns=100000;std::vector<Eigen::Triplet<double,std::int64_t>> entries;entries.reserve(3*ns);for(Eigen::Index i=0;i<ns;++i){entries.emplace_back(i,i,2.2);if(i){entries.emplace_back(i,i-1,-1.);entries.emplace_back(i-1,i,-1.);}}Sp S(ns,ns);S.setFromTriplets(entries.begin(),entries.end());S.makeCompressed();Eigen::VectorXd xt(ns);for(Eigen::Index i=0;i<ns;++i)xt[i]=std::sin(.001*i);Eigen::VectorXd b=S*xt;Eigen::SparseLU<Sp,Eigen::COLAMDOrdering<std::int64_t>> sparse;sparse.compute(S);check(sparse.info()==Eigen::Success,"sparse factor");Eigen::VectorXd xs=sparse.solve(b);double sparse_error=(xs-xt).norm()/xt.norm();check(sparse_error<1e-12,"sparse solve");std::cout<<"{\"case\":\"sparse_100k\",\"error\":"<<sparse_error<<",\"nnz\":"<<S.nonZeros()<<"}\n";
  Mat K=Mat::Identity(20,20)*2.;for(int i=1;i<20;++i)K(i,i-1)=K(i-1,i)=-.4;Mat C=Mat::Zero(3,20);C(0,0)=1.;C(1,7)=1.;C(2,19)=1.;Mat kkt=Mat::Zero(23,23);kkt.topLeftCorner(20,20)=K;kkt.topRightCorner(20,3)=C.transpose();kkt.bottomLeftCorner(3,20)=C;Eigen::VectorXd t=Eigen::VectorXd::LinSpaced(23,.1,2.);auto kk=checked_solve(kkt,kkt*t);check(kk.has_value()&&(kk.value()-t).norm()<1e-12,"indefinite KKT");
  Mat M=Mat::Identity(20,20);M.diagonal().array()+=Eigen::VectorXd::LinSpaced(20,0.,.5).array();K(0,0)=-1.;Eigen::GeneralizedSelfAdjointEigenSolver<Mat> modes(K,M);check(modes.info()==Eigen::Success&&modes.eigenvalues()[0]<0,"negative stiffness preserved");double modalres=(K*modes.eigenvectors()-M*modes.eigenvectors()*modes.eigenvalues().asDiagonal()).norm()/K.norm();check(modalres<1e-12,"generalized original residual");
  Mat P=Mat::Zero(64,64);for(int j=0;j<32;++j){double rho=1.01-.01*j,ang=.1+.02*j;P.block<2,2>(2*j,2*j)<<rho*std::cos(ang),-rho*std::sin(ang),rho*std::sin(ang),rho*std::cos(ang);if(j<31)P.block<2,2>(2*j,2*j+2).setConstant(.02);}test_arnoldi(P,6,24);
  CM Z=CM::Zero(32,32);for(int i=0;i<32;++i){Z(i,i)={1.05-.02*i,.05+.01*i};if(i<31)Z(i,i+1)={.02,-.01};}Spectra::DenseGenMatProd<std::complex<double>> cop(Z);Spectra::GenEigsSolver<decltype(cop)> ce(cop,4,16);ce.init();auto nc=ce.compute(Spectra::SortRule::LargestMagn,2000,1e-12);check(nc==4&&ce.info()==Spectra::CompInfo::Successful,"complex Arnoldi");double cr=0.;auto cv=ce.eigenvectors();auto cz=ce.eigenvalues();for(Eigen::Index j=0;j<cz.size();++j)cr=std::max(cr,(Z*cv.col(j)-cz[j]*cv.col(j)).norm()/(Z.norm()*cv.col(j).norm()));check(cr<1e-10,"complex original residual");std::cout<<"{\"case\":\"complex_Arnoldi_and_KKT\",\"complex_residual\":"<<cr<<",\"generalized_residual\":"<<modalres<<"}\n";
  std::cout<<"{\"passed\":true,\"eigen\":\""<<EIGEN_MAJOR_VERSION<<'.'<<EIGEN_MINOR_VERSION<<'.'<<EIGEN_PATCH_VERSION<<"\",\"index_bytes\":"<<sizeof(Eigen::Index)<<"}\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
