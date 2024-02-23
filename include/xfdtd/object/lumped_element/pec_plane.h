#ifndef _XFDTD_LIB_PEC_PLANE_H_
#define _XFDTD_LIB_PEC_PLANE_H_

#include <memory>

#include "xfdtd/object/object.h"

namespace xfdtd {

class PecPlane : public Object {
 public:
  PecPlane(std::string name, std::unique_ptr<Cube> cube);

  PecPlane(const PecPlane &) = delete;

  PecPlane(PecPlane &&) noexcept = default;

  PecPlane &operator=(const PecPlane &) = delete;

  PecPlane &operator=(PecPlane &&) noexcept = default;

  ~PecPlane() override = default;

  private:

};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_PEC_PLANE_H_
