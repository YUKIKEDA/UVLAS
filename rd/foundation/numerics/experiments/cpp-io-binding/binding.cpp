#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <nlohmann/json.hpp>
#include <vector>
#include <memory>
#include <set>
#include <cmath>
namespace py=pybind11;
struct Result{std::vector<double> values;};
std::string validate(const std::string& text){
  // Reject duplicates during SAX-style callbacks, before DOM overwrites can hide them.
  std::vector<std::set<std::string>> keys;
  auto callback=[&](int,nlohmann::json::parse_event_t event,nlohmann::json& item){
    using E=nlohmann::json::parse_event_t;
    if(event==E::object_start)keys.emplace_back();
    if(event==E::key&&!keys.back().insert(item.get<std::string>()).second)throw std::invalid_argument("InvalidInput:duplicate_key");
    if(event==E::object_end)keys.pop_back();return true;
  };
  auto j=nlohmann::json::parse(text,callback);if(!j.is_object())throw std::invalid_argument("InvalidInput:object");
  for(auto it=j.begin();it!=j.end();++it)if(it.key()!="schema"&&it.key()!="id"&&it.key()!="scale")throw std::invalid_argument("InvalidInput:unknown_key");
  if(j.at("schema")!="uvlas.design-io/1"||!j.at("id").is_string())throw std::invalid_argument("InvalidInput:schema_or_id");
  auto s=j.at("id").get<std::string>();if(s.empty()||s.front()=='0'||s.find_first_not_of("0123456789")!=std::string::npos)throw std::invalid_argument("InvalidInput:id");
  std::size_t consumed=0;auto id=std::stoull(s,&consumed);if(!id||consumed!=s.size())throw std::invalid_argument("InvalidInput:id");
  double scale=j.at("scale").get<double>();if(!std::isfinite(scale)||scale<=0.)throw std::invalid_argument("InvalidInput:scale");
  return j.dump();
}
PYBIND11_MODULE(_uvlas_design_binding,m){
  m.def("validate",&validate);
  py::class_<Result,std::shared_ptr<Result>>(m,"Result")
    .def("view",[](std::shared_ptr<Result> self){py::array_t<double> out({py::ssize_t(self->values.size())},{py::ssize_t(sizeof(double))},self->values.data(),py::cast(self));out.attr("setflags")(false);return out;});
  m.def("compute",[](std::int64_t n){if(n<0||n>1000000)throw std::invalid_argument("ResourceExceeded:n");auto result=std::make_shared<Result>();{py::gil_scoped_release release;result->values.resize(static_cast<std::size_t>(n));for(std::int64_t i=0;i<n;++i)result->values[static_cast<std::size_t>(i)]=double(i)*.25;}return result;});
}
