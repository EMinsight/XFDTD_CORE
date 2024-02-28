#ifndef _XFDTD_LIB_CORRECTOR_H_
#define _XFDTD_LIB_CORRECTOR_H_

namespace xfdtd {

class Corrector {
 public:
  virtual ~Corrector() = default;

  virtual void correctE() = 0;

  virtual void correctH() = 0;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_CORRECTOR_H_
