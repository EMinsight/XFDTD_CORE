#ifndef _XFDTD_CORE_EXCEPTION_H_
#define _XFDTD_CORE_EXCEPTION_H_

#include <exception>
#include <string>
#include <utility>

namespace xfdtd {
class XFDTDException : public std::exception {
 public:
  XFDTDException() = default;

  explicit XFDTDException(std::string message) : _message{std::move(message)} {}

  const char* what() const noexcept override { return _message.c_str(); }

 private:
  std::string _message{"XFDTD Exception"};
};
}  // namespace xfdtd

#endif  // _XFDTD_CORE_EXCEPTION_H_
