#ifndef _PIC_ACTION_RESERVE_HPP_
#define _PIC_ACTION_RESERVE_HPP_

#include "apt/action_base.hpp"
#include "field/field.hpp"

namespace pic {
  template < int D, typename R >
  struct ActionReserve {
  private:
    using Action_t = apt::ActionBase<D>;
    using field_type = field::Field<R,3,D>;

    std::vector<apt::array<apt::Range,D>> _range;
    std::vector<field_type> _Ersv;
    std::vector<field_type> _Brsv;

    std::vector<const Action_t*> _actions;
    const Action_t* _action_in_use = nullptr;

    inline bool _is_empty ( const apt::Index<D>& Ib, const apt::Index<D>& Ie ) const {
      for ( int i = 0; i < D; ++i ) {
        if ( Ib[i] >= Ie[i] ) return true;
      }
      return false;
    }

    std::optional<int> _is_nonempty_subset( const apt::Index<D>& Ib, const apt::Index<D>& Ie ) const {
      if ( _is_empty(Ib,Ie) ) return {};

      for ( int s = 0; s < _range.size(); ++s ) {
        bool is_sub = true;
        for ( int i = 0; i < D; ++i ) {
          if (  _range[s][i].begin() <= Ib[i] && Ie[i] <= _range[s][i].end() ) continue;
          else is_sub = false; break;
        }
        if ( is_sub ) return {s};
      }
      return {};
    }

    struct margin_iterator {
    private:
      int _index = 0;
      constexpr static int _end = D * 2; // each dimension has LFT and RGT
      const apt::ActionBase<D>& _act;

      int find_1st_nonempty_subset( int i_start ) const {
        while ( i_start < _end ) {
          int d = i_start / 2;
          int lr = i_start % 2;
          if ( _act[d].margin()[lr] < 1 ) ++i_start;
          else break;
        }
        return i_start;
      }

    public:
      using difference_type = void;
      using value_type = apt::pair<apt::Index<D>>;
      using reference = value_type;
      using pointer = void;
      using iterator_category = std::forward_iterator_tag;

      constexpr margin_iterator(const apt::ActionBase<D>& act) : _act(act) {}

      constexpr bool operator!=( int end ) const noexcept {
        return _index != end;
      }

      constexpr margin_iterator& operator++() noexcept {
        _index = find_1st_nonempty_subset(_index + 1);
        return *this;
      }

      constexpr margin_iterator operator++(int) noexcept {
        auto res = *this; ++(*this); return res;
      }

      constexpr reference operator*() noexcept {
        apt::Index<D> b, e;
        int d = _index / 2;
        int lr = _index % 2;
        b[d] = ( LFT == lr ) ? _act[d].far_begin() : _act[d].end();
        e[d] = b[d] + _act[d].margin()[lr]; // POLEDANCE check the logic
        for ( int i = 0; i < D; ++i ) { // set other directions
          if ( i == d ) continue;
          // if i is lower than d, use bulk + guard. This ensures no overlapping
          b[i] = _act[i].begin() - (i < d) * _act[i].margin()[LFT];
          e[i] = _act[i].end() + (i < d) * _act[i].margin()[RGT];
        }
        return {b,e};
      }

      constexpr margin_iterator begin() const noexcept {
        margin_iterator res = *this;
        // NOTE this line is essential
        res._index = find_1st_nonempty_subset(res._index);
        return res;
      }
      constexpr int end() const noexcept { return _end; }

    };

  public:
    void init (const std::vector<const Action_t*>& actions) { _actions = actions; }

    void reserve( const field_type& E, const field_type& B ) {
      for ( const auto* p : _actions ) {
        if ( !p ) continue;
        for ( auto [b,e] : margin_iterator(*p) ) {
          // TODO optimize: only save the non-subset part
          if (  _is_empty(b,e) || _is_nonempty_subset(b,e) ) continue;
          // add new ranges to reserve
          _range.push_back(apt::make_range(b,e,0));
          _Ersv.push_back({_range.back()});
          _Brsv.push_back({_range.back()});
          for ( const auto& I : apt::Block(b,e) ) {
            for ( int C = 0; C < 3; ++ C ) {
              _Ersv.back()[C](I) = E[C](I);
              _Brsv.back()[C](I) = B[C](I);
            }
          }
        }
      }
    }

    void revert_to_prior( field_type& E, field_type& B, const Action_t& action ) {
      // revert guard cells to previous values
      _action_in_use = &action;

      for ( auto [b,e] : margin_iterator(action) ) {
        auto sub = _is_nonempty_subset(b,e);
        if ( sub ) {
          const int s = *sub;
          for ( const auto& I : apt::Block(b,e) ) {
            for ( int C = 0; C < 3; ++ C ) {
              std::swap(_Ersv[s][C](I), E[C](I) );
              std::swap(_Brsv[s][C](I), B[C](I) );
            }
          }
        }

      }

    }

    void back_to_current( field_type& E, field_type& B ) {
      if ( _action_in_use ) revert_to_prior(E,B,*_action_in_use);
      _action_in_use = nullptr;
    }
  };
}

#endif
