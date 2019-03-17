#ifndef _APT_PRIORITY_QUEUE_HPP_
#define _APT_PRIORITY_QUEUE_HPP_
#include <queue>

namespace apt {
  // provides indexing. Indexing need not be distinct
  template < typename T, class Compare = std::less<T>, typename Ind = int >
  class priority_queue {
  private:

    struct node_t {
      Ind i{};
      T val{};
    };

    struct node_compare {
    private:
      Compare _comp;
    public:
      constexpr bool operator()( const node_t& a, const node_t& b ) noexcept {
        return _comp(a.val, b.val) ? true : a.i < b.i; //NOTE index breaks ties
      }
    };

    std::priority_queue< node_t, std::vector<node_t>, node_compare > _pq{};

  public:
    auto top() const {
      node_t tp = _pq.top();
      return std::tuple<Ind, const T&>(tp.i, tp.val);
    }

    inline void pop() { _pq.pop(); }

    inline void push(Ind i, T val) { _pq.push({i,val});}

    inline bool empty() const { return _pq.empty(); }

  };
}

#endif
