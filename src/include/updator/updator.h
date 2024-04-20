#ifndef _XFDTD_CORE_UPDATOR_H_
#define _XFDTD_CORE_UPDATOR_H_

#include <xfdtd/calculation_param/calculation_param.h>
#include <xfdtd/common/index_task.h>
#include <xfdtd/electromagnetic_field/electromagnetic_field.h>
#include <xfdtd/grid_space/grid_space.h>

#include <memory>

namespace xfdtd {

class XFDTDUpdatorException : public XFDTDException {
 public:
  explicit XFDTDUpdatorException(const std::string& message)
      : XFDTDException(message) {}
};

class Updator {
 public:
  Updator(std::shared_ptr<const GridSpace> grid_space,
          std::shared_ptr<const CalculationParam> calculation_param,
          std::shared_ptr<EMF> emf, IndexTask task);

  Updator(const Updator&) = default;

  Updator(Updator&&) noexcept = default;

  Updator& operator=(const Updator&) = default;

  Updator& operator=(Updator&&) noexcept = default;

  virtual ~Updator() = default;

  virtual void updateE() = 0;

  virtual void updateH() = 0;

  auto task() const { return _task; }

  virtual std::string toString() const;

  // only work for shared memory model
  bool containXNEdge() const { return task().xRange().start() == 0; }

  // only work for shared memory model
  bool containYNEdge() const { return task().yRange().start() == 0; }

  // only work for shared memory model
  bool containZNEdge() const { return task().zRange().start() == 0; }

 protected:
  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<const CalculationParam> _calculation_param;
  std::shared_ptr<EMF> _emf;

  virtual void updateEEdge() = 0;

  virtual void updateHEdge() = 0;

 private:
  IndexTask _task;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_UPDATOR_H_
