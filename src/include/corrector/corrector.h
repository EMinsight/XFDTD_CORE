#ifndef _XFDTD_CORE_CORRECTOR_H_
#define _XFDTD_CORE_CORRECTOR_H_

#include <string>
namespace xfdtd {

class Corrector {
 public:
  virtual ~Corrector() = default;

  virtual void correctE() = 0;

  virtual void correctH() = 0;

  virtual std::string toString() const { return "Corrector"; }
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_CORRECTOR_H_
