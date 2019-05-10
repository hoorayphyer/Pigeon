#ifndef _APT_PRIORITY_QUEUE_HPP_
#define _APT_PRIORITY_QUEUE_HPP_
#include <queue>

namespace apt {
  // provides indexing. Indexing need not be distinct
  template < typename T, class Compare = std::less<T>, typename Ind = int >
  class priority_queue {
  public:
    struct node_t {
      Ind i{};
      T val{};
    };

  private:
    struct node_compare {
    private:
      Compare _comp;
    public:
      constexpr bool operator()( const node_t& a, const node_t& b ) noexcept {
        //NOTE index breaks ties. NOTE use >= so that the smaller index gets chosen first
        return a.val == b.val ? a.i >= b.i : _comp(a.val, b.val);
      }
    };

    std::priority_queue< node_t, std::vector<node_t>, node_compare > _pq{};

  public:
    decltype(auto) top() const { return _pq.top(); }

    inline void pop() { _pq.pop(); }

    inline void push(Ind i, T val) { _pq.push({i,val});}

    inline bool empty() const { return _pq.empty(); }

    inline auto size() const { return _pq.size(); }
  };
}

#endif
