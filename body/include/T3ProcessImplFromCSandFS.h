#pragma once
#ifndef T3PROCESSIMPLFROMFSANDCS_H
#define T3PROCESSIMPLFROMFSANDCS_H

namespace t3 {

template <typename CSImpl, typename FSImpl> class ProcessImplFromCSandFS {
public:
  using CS_t = CSImpl;
  using FS_t = FSImpl;
  ProcessImplFromCSandFS(CSImpl cs, FSImpl fs):fCS(CSImpl(cs)),fFS(FSImpl(fs)) {}
  template <typename Floating, typename... Args> auto GetCS(Floating e, Args&... args) const;
  template <typename... Args> auto GetFS(Args&... args) const;
  
protected:
  CSImpl fCS;
  FSImpl fFS;
};

template <typename CSImpl, typename FSImpl>
template <typename Floating, typename... Args>
auto ProcessImplFromCSandFS<CSImpl, FSImpl>::GetCS(Floating e, Args&... args) const {
  return fCS.GetCS(e, args...);
}

template <typename CSImpl, typename FSImpl>
template <typename... Args>
auto ProcessImplFromCSandFS<CSImpl, FSImpl>::GetFS(Args&... args) const {
  return fFS.GetFS(args...);
}

} // namespace t3

#endif // T3PROCESSIMPLFROMFSANDCS_H
