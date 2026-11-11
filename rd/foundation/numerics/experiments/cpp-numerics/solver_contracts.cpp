#include <cstdint>
#include <Eigen/Dense>
#include <unsupported/Eigen/IterativeSolvers>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <iostream>
#include <cmath>
#include <stdexcept>
void check(bool ok,const char* s){if(!ok)throw std::runtime_error(s);}
int main(){try{
  using Mat=Eigen::MatrixXd;using Vec=Eigen::VectorXd;using Big=boost::multiprecision::cpp_bin_float_100;
  Mat A=Mat::Identity(96,96)*3.;for(int i=0;i<96;++i)for(int j=0;j<96;++j)A(i,j)+=.05*std::sin(.03*i+.07*j)/(1.+std::abs(i-j));Vec truth=Vec::LinSpaced(96,.1,2.),b=A*truth;
  Eigen::GMRES<Mat,Eigen::DiagonalPreconditioner<double>> gmres;gmres.set_restart(64);gmres.setMaxIterations(512);gmres.setTolerance(1e-13);gmres.compute(A);Vec x=gmres.solve(b);double gr=(A*x-b).norm()/b.norm();check(gmres.info()==Eigen::Success&&gr<1e-12,"fixed preconditioner GMRES original residual");
  // Nonlinear implicit chart C(r,d)=d+d^3-sin(r)=0; analytic J and J*.
  double r=.4,d=.2;for(int i=0;i<20;++i)d-=(d+d*d*d-std::sin(r))/(1.+3*d*d);
  double jd=std::cos(r)/(1.+3*d*d),h=1e-5;
  auto dependent=[](double z){double y=.2;for(int i=0;i<20;++i)y-=(y+y*y*y-std::sin(z))/(1.+3*y*y);return y;};
  double fd=(dependent(r+h)-dependent(r-h))/(2*h);check(std::abs(fd-jd)<1e-9,"implicit chart J");
  double tangent=.7,w0=.2,w1=-.3;double lhs=tangent*w0+jd*tangent*w1,rhs=tangent*(w0+jd*w1);check(std::abs(lhs-rhs)<1e-16,"implicit chart adjoint");
  // Kept low word is essential for a small constrained difference and a nonlinear consumer.
  double hi=1.,lo=std::ldexp(1.,-54);Big precise=Big(hi)+Big(lo);Big exact_square=precise*precise-Big(1);
  double p=hi*hi,err=std::fma(hi,hi,-p);Big consumed=Big(p-1.)+Big(err)+Big(2)*Big(hi)*Big(lo)+Big(lo)*Big(lo);
  check(abs(consumed-exact_square)<Big("1e-90"),"low word nonlinear consumption");
  double projected=hi+lo;check(projected==1.&&Big(projected)-1==0&&precise-1==Big(lo),"projection loses representable contract");
  bool accepted_projection=abs((Big(projected)-1)-Big(lo))/Big(lo)<Big("1e-11");check(!accepted_projection,"must reject uncertified projection");
  std::cout.precision(17);std::cout<<"{\"passed\":true,\"gmres_original_residual\":"<<gr<<",\"gmres_iterations\":"<<gmres.iterations()<<",\"implicit_chart_derivative_difference\":"<<std::abs(fd-jd)<<",\"kept_low_word_nonlinear_consumer\":true,\"uncertified_projection_rejected\":true}\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
