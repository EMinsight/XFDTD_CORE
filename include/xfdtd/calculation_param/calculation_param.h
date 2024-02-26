#ifndef _XFDTD_LIB_CALCULATION_PARAM_H_
#define _XFDTD_LIB_CALCULATION_PARAM_H_

#include <xfdtd/exception/exception.h>
#include <xfdtd/grid_space/grid_space.h>

#include <memory>
#include <vector>

#include "xfdtd/coordinate_system/coordinate_system.h"

namespace xfdtd {

class TimeParam {
 public:
  TimeParam(double dt, std::size_t size, std::size_t start_time_step = 0);

  TimeParam(const TimeParam&) = default;

  TimeParam(TimeParam&&) noexcept = default;

  TimeParam& operator=(const TimeParam&) = default;

  TimeParam& operator=(TimeParam&&) noexcept = default;

  virtual ~TimeParam() = default;

  double dt() const;

  std::size_t startTimeStep() const;

  std::size_t size() const;

  std::size_t endTimeStep() const;

  std::size_t currentTimeStep() const;

  std::size_t remainingTimeStep() const;

  virtual void nextStep();

  virtual void reset();

  xt::xarray<double> eTime() const;

  xt::xarray<double> hTime() const;

 private:
  double _dt;
  std::size_t _start_time_step;
  std::size_t _size;
  std::size_t _current_time_step;
};

// forward declare
class Material;
class MaterialParam {
 public:
  enum class Attribute { EPSILON, MU, SIGMA_E, SIGMA_M };

  enum class Property {
    EPS_X,
    EPS_Y,
    EPS_Z,
    MU_X,
    MU_Y,
    MU_Z,
    SIGMA_E_X,
    SIGMA_E_Y,
    SIGMA_E_Z,
    SIGMA_M_X,
    SIGMA_M_Y,
    SIGMA_M_Z
  };

  static Attribute fromPropertyToAttribute(Property property);

  static Axis::XYZ fromPropertyToXYZ(Property property);

  static Property fromAttributeAndXYZ(Attribute attribute, Axis::XYZ xyz);

  MaterialParam() = default;

  MaterialParam(const MaterialParam&) = default;

  MaterialParam(MaterialParam&&) noexcept = default;

  MaterialParam& operator=(const MaterialParam&) = default;

  MaterialParam& operator=(MaterialParam&&) noexcept = default;

  ~MaterialParam() = default;

  const xt::xarray<double>& epsX() const;

  const xt::xarray<double>& epsY() const;

  const xt::xarray<double>& epsZ() const;

  const xt::xarray<double>& muX() const;

  const xt::xarray<double>& muY() const;

  const xt::xarray<double>& muZ() const;

  const xt::xarray<double>& sigmaEX() const;

  const xt::xarray<double>& sigmaEY() const;

  const xt::xarray<double>& sigmaEZ() const;

  const xt::xarray<double>& sigmaMX() const;

  const xt::xarray<double>& sigmaMY() const;

  const xt::xarray<double>& sigmaMZ() const;

  const xt::xarray<double>& property(Property property) const;

  const xt::xarray<double>& property(Attribute attribute, Axis::XYZ xyz) const;

  const auto& materialArray() const;

  xt::xarray<double>& epsX();

  xt::xarray<double>& epsY();

  xt::xarray<double>& epsZ();

  xt::xarray<double>& muX();

  xt::xarray<double>& muY();

  xt::xarray<double>& muZ();

  xt::xarray<double>& sigmaEX();

  xt::xarray<double>& sigmaEY();

  xt::xarray<double>& sigmaEZ();

  xt::xarray<double>& sigmaMX();

  xt::xarray<double>& sigmaMY();

  xt::xarray<double>& sigmaMZ();

  xt::xarray<double>& property(Property property);

  xt::xarray<double>& property(Attribute attribute, Axis::XYZ xyz);

  auto&& materialArray();

  void addMaterial(std::shared_ptr<Material> material);

 private:
  xt::xarray<double> _eps_x;
  xt::xarray<double> _eps_y;
  xt::xarray<double> _eps_z;
  xt::xarray<double> _mu_x;
  xt::xarray<double> _mu_y;
  xt::xarray<double> _mu_z;
  xt::xarray<double> _sigma_e_x;
  xt::xarray<double> _sigma_e_y;
  xt::xarray<double> _sigma_e_z;
  xt::xarray<double> _sigma_m_x;
  xt::xarray<double> _sigma_m_y;
  xt::xarray<double> _sigma_m_z;

  std::vector<std::shared_ptr<Material>> _materials;
};

inline const auto& MaterialParam::materialArray() const { return _materials; }

inline auto&& MaterialParam::materialArray() { return _materials; }

class FDTDUpdateCoefficient {
 public:
  FDTDUpdateCoefficient() = default;

  FDTDUpdateCoefficient(const FDTDUpdateCoefficient&) = default;

  FDTDUpdateCoefficient(FDTDUpdateCoefficient&&) noexcept = default;

  FDTDUpdateCoefficient& operator=(const FDTDUpdateCoefficient&) = default;

  FDTDUpdateCoefficient& operator=(FDTDUpdateCoefficient&&) noexcept = default;

  ~FDTDUpdateCoefficient() = default;

  const xt::xarray<double>& cexe() const;

  const xt::xarray<double>& cexhy() const;

  const xt::xarray<double>& cexhz() const;

  const xt::xarray<double>& ceye() const;

  const xt::xarray<double>& ceyhz() const;

  const xt::xarray<double>& ceyhx() const;

  const xt::xarray<double>& ceze() const;

  const xt::xarray<double>& cezhx() const;

  const xt::xarray<double>& cezhy() const;

  const xt::xarray<double>& chxh() const;

  const xt::xarray<double>& chxey() const;

  const xt::xarray<double>& chxez() const;

  const xt::xarray<double>& chyh() const;

  const xt::xarray<double>& chyez() const;

  const xt::xarray<double>& chyex() const;

  const xt::xarray<double>& chzh() const;

  const xt::xarray<double>& chzex() const;

  const xt::xarray<double>& chzey() const;

  xt::xarray<double>& cexe();

  xt::xarray<double>& cexhy();

  xt::xarray<double>& cexhz();

  xt::xarray<double>& ceye();

  xt::xarray<double>& ceyhz();

  xt::xarray<double>& ceyhx();

  xt::xarray<double>& ceze();

  xt::xarray<double>& cezhx();

  xt::xarray<double>& cezhy();

  xt::xarray<double>& chxh();

  xt::xarray<double>& chxey();

  xt::xarray<double>& chxez();

  xt::xarray<double>& chyh();

  xt::xarray<double>& chyez();

  xt::xarray<double>& chyex();

  xt::xarray<double>& chzh();

  xt::xarray<double>& chzex();

  xt::xarray<double>& chzey();

 private:
  xt::xarray<double> _cexe;
  xt::xarray<double> _cexhy;
  xt::xarray<double> _cexhz;
  xt::xarray<double> _ceye;
  xt::xarray<double> _ceyhz;
  xt::xarray<double> _ceyhx;
  xt::xarray<double> _ceze;
  xt::xarray<double> _cezhx;
  xt::xarray<double> _cezhy;

  xt::xarray<double> _chxh;
  xt::xarray<double> _chxey;
  xt::xarray<double> _chxez;
  xt::xarray<double> _chyh;
  xt::xarray<double> _chyez;
  xt::xarray<double> _chyex;
  xt::xarray<double> _chzh;
  xt::xarray<double> _chzex;
  xt::xarray<double> _chzey;
};

class XFDTDCalculationParamException : public XFDTDException {
 public:
  explicit XFDTDCalculationParamException(std::string message);

  ~XFDTDCalculationParamException() override = default;
};

class CalculationParam {
 public:
  inline static constexpr double DEFAULT_CFL{0.98};

  CalculationParam(const GridSpace* grid_space, std::size_t start_time_step,
                   std::size_t time_size, double cfl = DEFAULT_CFL);

  ~CalculationParam() = default;

  double cfl() const;

  const std::unique_ptr<TimeParam>& timeParam() const;

  const std::unique_ptr<MaterialParam>& materialParam() const;

  const std::unique_ptr<FDTDUpdateCoefficient>& fdtdCoefficient() const;

  std::unique_ptr<TimeParam>& timeParam();

  std::unique_ptr<MaterialParam>& materialParam();

  std::unique_ptr<FDTDUpdateCoefficient>& fdtdCoefficient();

  void generateMaterialSpaceParam(const GridSpace* grid_space);

  void calculateCoefficient(const GridSpace* grid_space);

  static double calculateDt(double cfl, double dx, double dy, double dz);

 private:
  double _cfl;
  std::unique_ptr<TimeParam> _time_param;
  std::unique_ptr<MaterialParam> _material_param;
  std::unique_ptr<FDTDUpdateCoefficient> _fdtd_coefficient;

  void allocate(const GridSpace* grid_space);
};

}  // namespace xfdtd

#endif  // _XFDTD_LIB_CALCULATION_PARAM_H_
