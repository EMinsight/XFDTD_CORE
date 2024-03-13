#ifndef _XFDTD_CORE_RANDOM_H_
#define _XFDTD_CORE_RANDOM_H_

#include <random>
#include <string>

namespace xfdtd {

inline std::string randomString(size_t length) {
  static const std::string characters =
      "0123456789"
      "abcdefghijklmnopqrstuvwxyz"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  static std::mt19937 rng{std::random_device{}()};
  static std::uniform_int_distribution<std::string::size_type> dist(
      0, characters.length() - 1);

  std::string result;
  result.reserve(length);
  for (size_t i = 0; i < length; ++i) {
    result += characters[dist(rng)];
  }
  return result;
}

}  // namespace xfdtd

#endif  // _XFDTD_CORE_RANDOM_H_
