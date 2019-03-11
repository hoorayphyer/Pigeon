#ifndef  _APT_NUMERIC_HPP_
#define  _APT_NUMERIC_HPP_

#include "apt/vec_expression.hpp"
#include "apt/type_traits.hpp"

namespace apt {

  enum class BinaryOps : int {ADD=0, SUB, MUL, DIV};

  template < BinaryOps Op, typename E1, typename E2,
             typename T = most_precise_t< element_t<E1>, element_t<E2> > >
  struct FCompWise_Vec_Vec : public VecExpression<FCompWise_Vec_Vec<Op,E1,E2,T>, T, false> {
  private:
    const E1& _e1;
    const E2& _e2;
  public:
    using element_type = T;
    // NOTE use the smaller size here
    static constexpr int NDim = (E1::NDim < E2::NDim) ? E1::NDim : E2::NDim;

    constexpr FCompWise_Vec_Vec( const E1& e1, const E2& e2 ) noexcept
      : _e1(e1), _e2(e2){}

    template < int I >
    constexpr T v() const noexcept {
      static_assert(I < NDim);
      if constexpr ( Op == BinaryOps::ADD ) return _e1.template v<I>() + _e2.template v<I>();
      else if ( Op == BinaryOps::SUB ) return _e1.template v<I>() - _e2.template v<I>();
      else if ( Op == BinaryOps::MUL ) return _e1.template v<I>() * _e2.template v<I>();
      else if ( Op == BinaryOps::DIV ) return _e1.template v<I>() / _e2.template v<I>();
    }
  };

  template < BinaryOps Op, typename E, typename Real,
             typename T = most_precise_t< element_t<E>, Real > >
  struct FCompWise_Vec_Sca : public VecExpression<FCompWise_Vec_Sca<Op,E,Real,T>, T, false> {
  private:
    const E& _e;
    const Real& _t;

  public:
    using element_type = T;
    static constexpr int NDim = E::NDim;

    constexpr FCompWise_Vec_Sca( const E& e, const Real& t ) noexcept
      : _e(e), _t(t) {}

    template < int I >
    constexpr T v() const noexcept {
      if constexpr ( Op == BinaryOps::ADD ) return _e.template v<I>() + _t;
      else if ( Op == BinaryOps::SUB ) return _e.template v<I>() - _t;
      else if ( Op == BinaryOps::MUL ) return _e.template v<I>() * _t;
      else if ( Op == BinaryOps::DIV ) return _e.template v<I>() / _t;
    }
  };

}

namespace {
#define vec_def_op(_OP_, _OP_NAME_)                                     \
  template <typename E1, typename E2>                                   \
  constexpr auto operator _OP_ (const apt::VecExpression<E1>& e1,       \
                              const apt::VecExpression<E2>& e2)         \
    noexcept {                                                          \
    return                                                              \
      apt::FCompWise_Vec_Vec< apt::BinaryOps:: _OP_NAME_,               \
                              apt::VecExpression<E1>,                   \
                              apt::VecExpression<E2> > (e1, e2);        \
  }                                                                     \
                                                                        \
  template <typename E, typename Real>                                  \
  constexpr std::enable_if_t< std::is_arithmetic_v<Real>,               \
                              apt::FCompWise_Vec_Sca                    \
                              <apt::BinaryOps::_OP_NAME_,               \
                               apt::VecExpression<E>,                   \
                               Real > >                                 \
  operator _OP_ (const apt::VecExpression<E>& e, const Real& t)         \
    noexcept {                                                          \
    return                                                              \
      apt::FCompWise_Vec_Sca<apt::BinaryOps::_OP_NAME_,                 \
                             apt::VecExpression<E>,                     \
                             Real >(e, t);                              \
  }                                                                     \

  vec_def_op(+, ADD);
  vec_def_op(-, SUB);
  vec_def_op(*, MUL);
  vec_def_op(/, DIV);

  template <typename E1, typename E2 >
  constexpr E1& operator+= ( apt::VecExpression<E1>& v1, const apt::VecExpression<E2>& v2 ) noexcept {
    apt::foreach<0, E1::NDim>
      ( []( auto& x, const auto& y ) noexcept { x += y; }, v1, v2 );
    return static_cast<E1&>(v1);
  }

  template <typename E1, typename E2 >
  constexpr E1& operator-= ( apt::VecExpression<E1>& v1, const apt::VecExpression<E2>& v2 ) noexcept {
    apt::foreach<0, E1::NDim>
      ( []( auto& x, const auto& y ) noexcept { x -= y; }, v1, v2 );
    return static_cast<E1&>(v1);
  }

  template <typename E, typename Real = typename E::element_type>
  constexpr E& operator*= ( apt::VecExpression<E>& v, const Real& t ) noexcept {
    apt::foreach<0, E::NDim>
      ( [&t]( auto& x ) noexcept { x *= t; }, v );
    return static_cast<E&>(v);
  }

  template <typename E, typename Real = typename E::element_type>
  constexpr E& operator/= ( apt::VecExpression<E>& v, const Real& t ) noexcept {
    apt::foreach<0, E::NDim>
      ( [&t]( auto& x ) noexcept { x /= t; }, v );
    return static_cast<E&>(v);
  }
}

namespace apt {
  template <typename E1, typename E2, typename T = typename E1::element_type>
  struct VecCross : public VecExpression<VecCross<E1,E2,T>, T, false> {
  private:
    const E1& _e1;
    const E2& _e2;

  public:
    using element_type = T;
    static constexpr int NDim = E1::NDim;

    constexpr VecCross( const E1& e1, const E2& e2 ) noexcept
      : _e1(e1), _e2(e2) {}

    template < int I >
    constexpr T v() const noexcept {
      return _e1.template v<(I+1)%3>() * _e2.template v<(I+2)%3>()
        - _e1.template v<(I+2)%3>() * _e2.template v<(I+1)%3>();
    }
  };

  template < typename E1, typename E2 >
  constexpr auto cross( const E1& e1, const E2& e2 ) noexcept {
    return VecCross<E1,E2>(e1, e2);
  }
}

#include <cmath>
namespace apt {
  template < int I, typename E1, typename E2 >
  constexpr auto dot ( const VecExpression<E1>& v1, const VecExpression<E2>& v2 ) noexcept {
  }

  template < typename E1, typename E2, int I = 0 >
  constexpr auto dot ( const VecExpression<E1>& v1, const VecExpression<E2>& v2 ) noexcept {
    constexpr int END = E1::NDim;
    if constexpr ( I == END ) return 0;
    else return v1.template v<I>() * v2.template v<I>() + dot<E1, E2, I+1>(v1, v2);
  }

  template < typename E >
  constexpr auto sqabs( const VecExpression<E>& vec ) noexcept {
    return dot(vec, vec);
  }

  template < typename E >
  constexpr auto abs( const VecExpression<E>& vec ) noexcept {
    constexpr int N = E::NDim;
    if constexpr ( N == 1 ) return std::abs( vec.template v<0>() );
    else return std::sqrt( sqabs(vec) );
  }

}


#endif
