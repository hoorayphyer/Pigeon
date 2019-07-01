if (${ENABLE_TESTS})
  test_gen(manifold curvilinear) # DONE
  # test_gen(manifold grid) # TODO

  msperf_gen(manifold curvilinear) # TODO
endif()
