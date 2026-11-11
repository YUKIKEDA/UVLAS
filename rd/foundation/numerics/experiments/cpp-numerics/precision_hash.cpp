#include <cstdint>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/hash2/sha2.hpp>
#include <boost/hash2/digest.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <array>
void check(bool x,const char* message){if(!x)throw std::runtime_error(message);}
int main(){try{
  using boost::multiprecision::cpp_int;using boost::multiprecision::cpp_rational;using Big=boost::multiprecision::cpp_bin_float_100;
  double a=1.+std::ldexp(1.,-27),b=1.-std::ldexp(1.,-27);double product=a*b;double low=std::fma(a,b,-product);double compensated=(product-1.)+low;Big oracle=Big(a)*Big(b)-1;check(compensated==-std::ldexp(1.,-54)&&abs(Big(compensated)-oracle)<Big("1e-90"),"error-free product");
  double sum=0.,correction=0.;for(double value:std::array<double,3>{1.e16,1.,-1.e16}){double next=sum+value;if(std::abs(sum)>=std::abs(value))correction+=(sum-next)+value;else correction+=(value-next)+sum;sum=next;}check(sum+correction==1.,"compensated cancellation");
  cpp_rational tiny=cpp_rational(1)/cpp_rational(cpp_int(1)<<180);cpp_rational determinant=(1+tiny)-1;check(determinant!=0&&(1+tiny).convert_to<double>()==1.,"exact rank versus float collapse");
  cpp_int scale=cpp_int(1)<<256;cpp_int radicand=3*scale*scale/2;cpp_int root=sqrt(radicand);check(root*root<=radicand&&(root+1)*(root+1)>radicand,"integer square root");cpp_rational lo=cpp_rational(root)/cpp_rational(scale),hi=cpp_rational(root+1)/cpp_rational(scale),target=cpp_rational(3)/2;check(lo*lo<=target&&hi*hi>=target&&hi-lo==cpp_rational(1)/cpp_rational(scale),"rational outward sqrt enclosure");
  boost::hash2::sha2_256 hash;hash.update("abc",3);auto digest=boost::hash2::to_string(hash.result());check(digest=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad","SHA256 known answer");
  std::cout.precision(17);std::cout<<"{\"passed\":true,\"two_product_low\":"<<low<<",\"compensated_sum\":"<<sum+correction<<",\"exact_rank_preserved\":true,\"sqrt_enclosure_bits\":256,\"sha256_abc\":\""<<digest<<"\"}\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
