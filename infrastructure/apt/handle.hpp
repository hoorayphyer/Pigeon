#pragma once

#include <memory>

namespace apt {
template <typename RawHdl, void (*Deleter)(RawHdl*), RawHdl (*Default)()>
struct Handle {
 private:
  using shared_ptr = std::shared_ptr<RawHdl>;
  shared_ptr _ptr;
  RawHdl _default = Default();
  RawHdl _fallback = Default();

 public:
  Handle() = default;
  Handle(const Handle&) = default;
  Handle(Handle&& other) noexcept = default;
  virtual ~Handle() = default;

  Handle& operator=(const Handle&) noexcept = default;
  Handle& operator=(Handle&& other) noexcept = default;

  inline void reset() { _ptr.reset(); }

  inline void reset(RawHdl* raw_handle) { _ptr.reset(raw_handle, Deleter); }

  inline void set_fallback(RawHdl fallback) { _fallback = fallback; }

  inline void reset_fallback() { _fallback = _default; }

  operator RawHdl() const noexcept { return use_count() ? *_ptr : _fallback; }

  operator RawHdl*() const noexcept { return _ptr.get(); }

  operator bool() const noexcept { return static_cast<bool>(_ptr); }

  long use_count() const noexcept { return _ptr.use_count(); }
};
}  // namespace apt
