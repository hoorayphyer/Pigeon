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

    std::vector<apt::pair<apt::Index<D>>> _IbIe;
    std::vector<apt::array<std::vector<R>,3>> _Ersv;
    std::vector<apt::array<std::vector<R>,3>> _Brsv;

    std::vector<const Action_t*> _actions;
    const Action_t* _action_in_use = nullptr;

    std::optional<std::optional<int>> _is_subset( const apt::Index<D>& Ib, const apt::Index<D>& Ie ) const {
      // two layers of optional because there is not-a-subset, a-nonempty-subset, and empty subset
      for ( int s = 0; s < _IbIe.size(); ++s ) {
        const auto& [b,e] = _IbIe[s];

        bool is_sub = true;
        for ( int i = 0; i < D; ++i ) {
          // if Ib Ie forms an empty set, return empty optional
          if ( Ib[i] >= Ie[i] ) return {{}}; // TODO double check
          if (  b[i] <= Ib[i] && Ie[i] <= e[i] ) continue;
          is_sub = false; break;
        }
        if ( is_sub ) return {{s}};
      }
      return {};
    }

    struct guard_range_iterator {
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

      constexpr guard_range_iterator(const apt::ActionBase<D>& act) : _act(act) {}

      constexpr bool operator!=( int end ) const noexcept {
        return _index != end;
      }

      constexpr guard_range_iterator& operator++() noexcept {
        _index = find_1st_nonempty_subset(_index + 1);
        return *this;
      }

      constexpr guard_range_iterator operator++(int) noexcept {
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

      constexpr guard_range_iterator begin() const noexcept {
        guard_range_iterator res = *this;
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
        for ( auto [b,e] : guard_range_iterator(*p) ) {
          // TODO optimize: only save the non-subset part
          if ( _is_subset(b,e) ) continue;
          _IbIe.push_back({b,e});
          _Ersv.push_back({});
          _Brsv.push_back({});
          for ( const auto& I : apt::Block(b,e) ) {
            for ( int C = 0; C < 3; ++ C ) {
              _Ersv.back()[C].push_back(E[C](I));
              _Brsv.back()[C].push_back(B[C](I));
            }
          }
        }
      }
    }

    void revert_to_prior( field_type& E, field_type& B, const Action_t& action ) {
      // revert guard cells to previous values
      _action_in_use = &action;

      // FIXME
      for ( auto [b,e] : guard_range_iterator(action) ) {
        auto sub = _is_subset(b,e);
        if ( sub && !*sub ) { // non empty subset
          const int s = *(*sub);
          const auto& [Ib, Ie] = _IbIe[s];
          int shift = b[D-1] - Ib[D-1];
          for ( int i = D-2 ; i > -1; --i ) {
            shift = shift * ( Ie[i] - Ib[i] ) + b[i] - Ib[i];
          }
          for ( auto I : apt::Block({},e-b) ) {
            int i_rsv = I[D-1];
            for ( int i = D-2; i > -1; --i ) i_rsv = i_rsv * (Ie[i]-Ib[i]) + I[i];
            i_rsv += shift;

            I += b;
            for ( int C = 0; C < 3; ++ C ) {
              std::swap(_Ersv[s][C][i_rsv], E[C](I) );
              std::swap(_Brsv[s][C][i_rsv], B[C](I) );
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
