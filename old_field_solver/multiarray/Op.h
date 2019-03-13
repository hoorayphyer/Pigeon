#ifndef _OP_H_
#define _OP_H_

template <typename T>
struct Op_Assign {
  inline void operator()(T& dest, const T& value) const { dest = value; }
};

template <typename T>
struct Op_PlusAssign {
  inline void operator()(T& dest, const T& value) const { dest += value; }
};

template <typename T>
struct Op_MinusAssign {
  inline void operator()(T& dest, const T& value) const { dest -= value; }
};

template <typename T>
struct Op_MultAssign {
  inline void operator()(T& dest, const T& value) const { dest *= value; }
};

template <typename T>
struct Op_AssignConst {
  T _value;
  Op_AssignConst(const T& value) : _value(value) {}

  inline void operator()(T& dest) const { dest = _value; }
};

template <typename T>
struct Op_MultConst {
  T _value;
  Op_MultConst(const T& value) : _value(value) {}

  inline void operator()(T& dest) const { dest *= _value; }
};

template <typename T>
struct Op_PlusConst {
  T _value;
  Op_PlusConst(const T& value) : _value(value) {}

  inline void operator()(T& dest) const { dest += _value; }
};

template <typename T>
struct Op_MinusConst {
  T _value;
  Op_MinusConst(const T& value) : _value(value) {}

  inline void operator()(T& dest) const { dest -= _value; }
};

#endif  // --- ifndef _OP_H_
