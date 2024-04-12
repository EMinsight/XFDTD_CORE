#include <xfdtd/calculation_param/calculation_param.h>

#include <xtensor/xnpy.hpp>

namespace xfdtd {

auto FDTDUpdateCoefficient::save(const std::string& dir) const -> void {
  auto dir_path = std::filesystem::path(dir);

  if (!std::filesystem::exists(dir_path) ||
      !std::filesystem::is_directory(dir_path)) {
    std::filesystem::create_directories(dir_path);
  }

  xt::dump_npy(dir_path / "cexe.npy", _cexe);
  xt::dump_npy(dir_path / "cexhy.npy", _cexhy);
  xt::dump_npy(dir_path / "cexhz.npy", _cexhz);
  xt::dump_npy(dir_path / "ceye.npy", _ceye);
  xt::dump_npy(dir_path / "ceyhx.npy", _ceyhx);
  xt::dump_npy(dir_path / "ceyhz.npy", _ceyhz);
  xt::dump_npy(dir_path / "ceze.npy", _ceze);
  xt::dump_npy(dir_path / "cezhx.npy", _cezhx);
  xt::dump_npy(dir_path / "cezhy.npy", _cezhy);
  xt::dump_npy(dir_path / "chxh.npy", _chxh);
  xt::dump_npy(dir_path / "chxey.npy", _chxey);
  xt::dump_npy(dir_path / "chxez.npy", _chxez);
  xt::dump_npy(dir_path / "chyh.npy", _chyh);
  xt::dump_npy(dir_path / "chyez.npy", _chyez);
  xt::dump_npy(dir_path / "chyex.npy", _chyex);
  xt::dump_npy(dir_path / "chzh.npy", _chzh);
  xt::dump_npy(dir_path / "chzex.npy", _chzex);
  xt::dump_npy(dir_path / "chzey.npy", _chzey);
}

const Array3D<Real>& FDTDUpdateCoefficient::cexe() const { return _cexe; }

const Array3D<Real>& FDTDUpdateCoefficient::cexhy() const { return _cexhy; }

const Array3D<Real>& FDTDUpdateCoefficient::cexhz() const { return _cexhz; }

const Array3D<Real>& FDTDUpdateCoefficient::ceye() const { return _ceye; }

const Array3D<Real>& FDTDUpdateCoefficient::ceyhx() const { return _ceyhx; }

const Array3D<Real>& FDTDUpdateCoefficient::ceyhz() const { return _ceyhz; }

const Array3D<Real>& FDTDUpdateCoefficient::ceze() const { return _ceze; }

const Array3D<Real>& FDTDUpdateCoefficient::cezhx() const { return _cezhx; }

const Array3D<Real>& FDTDUpdateCoefficient::cezhy() const { return _cezhy; }

const Array3D<Real>& FDTDUpdateCoefficient::chxh() const { return _chxh; }

const Array3D<Real>& FDTDUpdateCoefficient::chxey() const { return _chxey; }

const Array3D<Real>& FDTDUpdateCoefficient::chxez() const { return _chxez; }

const Array3D<Real>& FDTDUpdateCoefficient::chyh() const { return _chyh; }

const Array3D<Real>& FDTDUpdateCoefficient::chyez() const { return _chyez; }

const Array3D<Real>& FDTDUpdateCoefficient::chyex() const { return _chyex; }

const Array3D<Real>& FDTDUpdateCoefficient::chzh() const { return _chzh; }

const Array3D<Real>& FDTDUpdateCoefficient::chzex() const { return _chzex; }

const Array3D<Real>& FDTDUpdateCoefficient::chzey() const { return _chzey; }

Array3D<Real>& FDTDUpdateCoefficient::cexe() { return _cexe; }

Array3D<Real>& FDTDUpdateCoefficient::cexhy() { return _cexhy; }

Array3D<Real>& FDTDUpdateCoefficient::cexhz() { return _cexhz; }

Array3D<Real>& FDTDUpdateCoefficient::ceye() { return _ceye; }

Array3D<Real>& FDTDUpdateCoefficient::ceyhx() { return _ceyhx; }

Array3D<Real>& FDTDUpdateCoefficient::ceyhz() { return _ceyhz; }

Array3D<Real>& FDTDUpdateCoefficient::ceze() { return _ceze; }

Array3D<Real>& FDTDUpdateCoefficient::cezhx() { return _cezhx; }

Array3D<Real>& FDTDUpdateCoefficient::cezhy() { return _cezhy; }

Array3D<Real>& FDTDUpdateCoefficient::chxh() { return _chxh; }

Array3D<Real>& FDTDUpdateCoefficient::chxey() { return _chxey; }

Array3D<Real>& FDTDUpdateCoefficient::chxez() { return _chxez; }

Array3D<Real>& FDTDUpdateCoefficient::chyh() { return _chyh; }

Array3D<Real>& FDTDUpdateCoefficient::chyez() { return _chyez; }

Array3D<Real>& FDTDUpdateCoefficient::chyex() { return _chyex; }

Array3D<Real>& FDTDUpdateCoefficient::chzh() { return _chzh; }

Array3D<Real>& FDTDUpdateCoefficient::chzex() { return _chzex; }

Array3D<Real>& FDTDUpdateCoefficient::chzey() { return _chzey; }

}  // namespace xfdtd
