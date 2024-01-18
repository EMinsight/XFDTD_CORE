#ifndef _XFDTD_LIB_CORRECTOR_H_
#define _XFDTD_LIB_CORRECTOR_H_

namespace xfdtd {

class Corrector {
 public:
  Corrector() = default;

  virtual ~Corrector() = default;

  virtual void correct() = 0;
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_CORRECTOR_H_
