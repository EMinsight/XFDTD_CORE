#include <xfdtd/monitor/movie_monitor.h>

#include <iomanip>
#include <memory>
#include <sstream>

namespace xfdtd {

MovieMonitor::MovieMonitor(std::unique_ptr<Monitor> frame,
                           std::size_t frame_interval, std::string name,
                           std::string output_dir_path)
    : Monitor{{nullptr}, std::move(name), std::move(output_dir_path)},
      _frame{std::move(frame)},
      _frame_interval{frame_interval},
      _frame_count{0} {}

void MovieMonitor::init(
    std::shared_ptr<const GridSpace> grid_space,
    std::shared_ptr<const CalculationParam> calculation_param,
    std::shared_ptr<const EMF> emf) {
  defaultInit(grid_space, calculation_param, emf);
  _frame->init(grid_space, calculation_param, emf);
  _frame->setOutputDir(outputDir());
}

void MovieMonitor::update() {
  if (_frame_count % _frame_interval == 0) {
    frame()->setName(formatFrameCount(_frame_count));
    frame()->update();
    frame()->output();
  }

  _frame_count++;
}

void MovieMonitor::output() {}

auto MovieMonitor::initParallelizedConfig() -> void {
  frame()->initParallelizedConfig();
}

std::size_t MovieMonitor::frameInterval() const { return _frame_interval; }

std::size_t MovieMonitor::frameCount() const { return _frame_count; }

const std::unique_ptr<Monitor>& MovieMonitor::frame() const { return _frame; }

std::unique_ptr<Monitor>& MovieMonitor::frame() { return _frame; }

void MovieMonitor::setFrameInterval(std::size_t frame_interval) {
  _frame_interval = frame_interval;
}

void MovieMonitor::setFrameCount(std::size_t frame_count) {
  _frame_count = frame_count;
}

std::string MovieMonitor::formatFrameCount(std::size_t frame_count) const {
  std::string frame_count_str{std::to_string(frame_count)};
  std::stringstream ss;
  ss << std::setw(5) << std::setfill('0') << frame_count_str;
  return ss.str();
}

std::string MovieMonitor::toString() const {
  std::stringstream ss;
  ss << "Movie Monitor: " << name() << "\n";
  ss << " Frame Interval: " << frameInterval() << "\n";
  ss << " Frame Count: " << frameCount() << "\n";
  ss << " Frame monitor: \n";
  ss << "  " << frame()->toString();
  return ss.str();
}

auto MovieMonitor::valid() const -> bool { return frame()->valid(); }

}  // namespace xfdtd
