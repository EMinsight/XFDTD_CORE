#ifndef _XFDTD_CORE_MOVIE_MONITOR_H_
#define _XFDTD_CORE_MOVIE_MONITOR_H_

#include <xfdtd/monitor/monitor.h>

#include <cstddef>
#include <memory>

namespace xfdtd {

class MovieMonitor : public xfdtd::Monitor {
 public:
  MovieMonitor(std::unique_ptr<Monitor> frame, std::size_t frame_interval,
               std::string name = "movie_monitor",
               std::string output_dir_path = "xfdtd_output");

  MovieMonitor(const MovieMonitor&) = delete;

  MovieMonitor(MovieMonitor&&) noexcept = default;

  MovieMonitor& operator=(const MovieMonitor&) = delete;

  MovieMonitor& operator=(MovieMonitor&&) noexcept = default;

  ~MovieMonitor() override = default;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<const CalculationParam> calculation_param,
            std::shared_ptr<const EMF> emf) override;

  void update() override;

  void output() override;

  std::size_t frameInterval() const;

  std::size_t frameCount() const;

  const std::unique_ptr<xfdtd::Monitor>& frame() const;

  std::unique_ptr<xfdtd::Monitor>& frame();

  void setFrameInterval(std::size_t frame_interval);

  void setFrameCount(std::size_t frame_count);

 private:
  std::unique_ptr<xfdtd::Monitor> _frame;

  std::size_t _frame_interval;

  std::size_t _frame_count;

  std::string formatFrameCount(std::size_t frame_count) const;
};

}  // namespace xfdtd

#endif  // _XFDTD_CORE_MOVIE_MONITOR_H_
