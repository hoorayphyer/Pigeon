#ifndef _PARTICLE_STATE_CODEC_HPP_
#define _PARTICLE_STATE_CODEC_HPP_

#include <type_traits>
#include <bitset>

namespace particle {
  enum class species : char {
    electron, positron, ion, photon
  };

  template < typename T, std::size_t NBits >
  class species_manager
    : private {
  private:
  public:
    template< species s >
    void set();
    void species() const; // decode species
  };

  enum class flag : int {
    empty = 0, secondary, ignore_force,
    ignore_deposit, annihilate, ignore_em,
    delimiter, traced
  };

  template < typename T, std::size_t NBits >
  class flagger {
  private:
    // T& _state;

  public:
    // flagger( T& s ) noexcept : _state(s) {}

    template < flag fl, flag... rest >
    inline void set() const noexcept {
      if constexpr ( flag::traced == fl ) {
          if ( !is<flag::traced>() ) encode_trace_id();
        }
      static_cast<bitset_t>(_state).set( static_cast<int>(fl) );
      if constexpr ( sizeof...(rest) > 0 ) set<rest...>();
    }

    template < flag fl, flag... rest >
    inline void reset() const noexcept {
      static_cast<bitset_t>(_state).reset( static_cast<int>(fl) );
      if constexpr ( sizeof...(rest) > 0 ) reset<rest...>();
    }

    template < flag fl, flag... rest >
    inline void flip() const noexcept {
      if constexpr ( flag::traced == fl ) {
          if ( !is<flag::traced>() ) encode_trace_id();
        }
      static_cast<bitset_t>(_state).flip( static_cast<int>(fl) );
      if constexpr ( sizeof...(rest) > 0 ) flip<rest...>();
    }

    template < flag fl >
    inline bool is() const noexcept {
      return static_cast<const bitset_t>(_state)[static_cast<int>(fl)];
    }


  };

  template < typename T, std::size_t NB_Rank, std::size_t NB_Count >
  class tracker {
  private:
    unsigned int trace_counter = 0; // TODOL how to avoid overflow because only 32 bits are assigned
    static const int rank = 0;

  public:
    // TODO check
    // TODOL how to check overflow
    static void encode_trace_id() {
      auto& state_bits = static_cast<bitset_t>(_state);
      const auto& rank_bits = static_cast<std::bitset<RankBitSize>>(rank);
      const auto& count_bits = static_cast<std::bitset<8 * sizeof(encoded_bits_t)-FlagBitSize-RankBitSize>>(trace_counter);
      // TODO impl
      // copy(rank_bits, state_bits[FlagBitSize]);
      // copy(count_bits, state_bits[FlagBitSize+RankBitSize]);

      ++trace_counter;
    }

  };




  template < typename T, int NBits = 8 * sizeof( std::remove_reference_t<T> ),
             int NB_species = 2, int NB_flag = 16, int NB_rank = 16,
             int NB_count = NBits - NB_species - NB_flag - NB_rank >
  struct state_codec : private species_manager<T, NB_species>
                     , private flagger<T, NB_flag>
                     , private tracker<T, NB_rank, NB_count> {
  private:
    using bitset_t = std::bitset<NBits>;

  protected:
    T state; // NOTE use T here so when T is nonref, _state will hold the data


  public:
    explicit state_codec( T my_state ) noexcept : state(my_state) {}


  };
}

using particle::flag;
using particle::species;


#endif
