#ifndef _PARTICLE_STATE_CODEC_HPP_
#define _PARTICLE_STATE_CODEC_HPP_

#include "types.hpp"
#include <bitset>

namespace particle {
  enum class flag : int {
    empty = 0, secondary, ignore_force,
    ignore_deposit, annihilate, ignore_em,
    delimiter, traced
  };

  struct state {
  private:
    static constexpr int FlagBitSize = 16;
    static constexpr int RankBitSize = 16; // for tracing

    static unsigned int trace_counter = 0; // TODOL how to avoid overflow because only 32 bits are assigned
    static const int rank = 0;

    using bitset_t = std::bitset<8 * sizeof(encoded_bits_t) >;

    // TODO check
    // TODOL how to check overflow
    static void encode_trace_id( encoded_bits_t& state ) {
      auto& state_bits = static_cast<bitset_t>(state);
      const auto& rank_bits = static_cast<std::bitset<RankBitSize>>(rank);
      const auto& count_bits = static_cast<std::bitset<8 * sizeof(encoded_bits_t)-FlagBitSize-RankBitSize>>(trace_counter);
      copy(rank_bits, state_bits[FlagBitSize]);
      copy(count_bits, state_bits[FlagBitSize+RankBitSize]);

      ++trace_counter;
    }

  public:
    template < flag f, flag... rest >
    static inline void set( encoded_bits_t& state ) noexcept {
      if constexpr ( flag::traced == f ) {
        if ( !is<flag::traced>(state) ) encode_trace_id(state);
      }
      static_cast<bitset_t>(state).set( static_cast<int>(fl) );
      if constexpr ( sizeof...(rest) > 0 ) set<rest...>(state);
    }

    template < flag... FLAG >
    static inline void reset( encoded_bits_t& state, flag fl, FLAG... rest ) noexcept {
      static_cast<bitset_t>(state).reset( static_cast<int>(fl) );
      if constexpr ( sizeof...(rest) > 0 ) reset(state, rest...);
    }

    template < flag... FLAG >
    static inline void flip( encoded_bits_t& state, flag fl, FLAG... rest ) noexcept {
      if constexpr ( flag::traced == f ) {
        if ( !is<flag::traced>(state) ) encode_trace_id(state);
      }
      static_cast<bitset_t>(state).flip( static_cast<int>(fl) );
      if constexpr ( sizeof...(rest) > 0 ) flip(state, rest...);
    }

    template < flag FLAG >
    static inline bool is( const encoded_bits_t& state ) noexcept {
      return static_cast<bitset_t>(state)[static_cast<int>(FLAG)];
    }
  };
}


#endif
