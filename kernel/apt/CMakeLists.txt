if (${ENABLE_TESTS})
  test_gen(apt array)
  test_gen(apt pair)
  test_gen(apt "foreach")
  test_gen(apt block)
  test_gen(apt vec)
  # test_gen(apt virtual_vec) # TODO
  test_gen(apt numeric)
  test_gen(apt handle)
  test_gen(apt priority_queue)
  test_gen(apt csi)

  msperf_gen(apt block)
endif()
