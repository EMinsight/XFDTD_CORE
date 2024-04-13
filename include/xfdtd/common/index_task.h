#ifndef __XFDTD_CORE_INDEX_TASK_H__
#define __XFDTD_CORE_INDEX_TASK_H__

#include <xfdtd/common/type_define.h>

#include <optional>

namespace xfdtd {

template <typename T>
class Range {
 public:
  Range() = default;

  Range(T start, T end) : _start(start), _end(end) {}

  auto operator+(T value) const { return Range{_start + value, _end + value}; }

  auto operator-(T value) const { return Range{_start - value, _end - value}; }

  auto start() const { return _start; }

  auto end() const { return _end; }

  auto size() const { return _end - _start; }

  auto valid() const { return _start < _end; }

  auto toString() const {
    std::stringstream ss;
    ss << "[" << _start << " " << _end << ")";
    return ss.str();
  }

 private:
  T _start{};
  T _end{};
};

template <typename T>
struct Task {
  Range<T> _x_range;
  Range<T> _y_range;
  Range<T> _z_range;

  auto toString() const {
    std::stringstream ss;
    ss << "Task: ";
    ss << "x: " << _x_range.toString() << " ";
    ss << "y: " << _y_range.toString() << " ";
    ss << "z: " << _z_range.toString();
    return ss.str();
  }

  auto xRange() const { return _x_range; }

  auto yRange() const { return _y_range; }

  auto zRange() const { return _z_range; }

  auto valid() const {
    return _x_range.valid() && _y_range.valid() && _z_range.valid();
  }
};

template <typename T>
inline auto makeRange(T start, T end) {
  return Range<T>{start, end};
}

template <typename T>
inline auto makeTask(const Range<T>& x_range, const Range<T>& y_range,
                     const Range<T>& z_range) {
  return Task<T>{x_range, y_range, z_range};
}

using IndexRange = Range<Index>;

using IndexTask = Task<Index>;

inline auto makeIndexRange(Index start, Index end) {
  return makeRange(start, end);
}

inline auto makeIndexTask(const Range<Index>& x_range,
                          const Range<Index>& y_range,
                          const Range<Index>& z_range) {
  return makeTask(x_range, y_range, z_range);
}

template <typename T>
inline bool intersected(const Range<T>& domain, const Range<T>& range) {
  return !(range.end() <= domain.start() || domain.end() <= range.start());
}

template <typename T>
inline bool intersected(const Task<T>& task_1, const Task<T>& task_2) {
  return intersected(task_1._x_range, task_2._x_range) &&
         intersected(task_1._y_range, task_2._y_range) &&
         intersected(task_1._z_range, task_2._z_range);
}

template <typename T>
inline std::optional<Task<T>> taskIntersection(const Task<T>& task1,
                                               const Task<T>& task2) {
  if (!intersected(task1, task2)) {
    return {};
  }

  auto x_range =
      Range<T>{std::max(task1.xRange().start(), task2.xRange().start()),
               std::min(task1.xRange().end(), task2.xRange().end())};
  auto y_range =
      Range<T>{std::max(task1.yRange().start(), task2.yRange().start()),
               std::min(task1.yRange().end(), task2.yRange().end())};
  auto z_range =
      Range<T>{std::max(task1.zRange().start(), task2.zRange().start()),
               std::min(task1.zRange().end(), task2.zRange().end())};
  return Task<T>{x_range, y_range, z_range};
}

}  // namespace xfdtd

#endif  // __XFDTD_CORE_INDEX_TASK_H__
