if (${ENABLE_TESTS})
  test_gen(apt array)
  test_gen(apt pair)
  test_gen(apt "foreach")
  test_gen(apt index)
  test_gen(apt block)
  test_gen(apt grid)
  test_gen(apt vec)
  test_gen(apt virtual_vec)
  test_gen(apt numeric)
  test_gen(apt handle)
  test_gen(apt priority_queue)
  test_gen(apt csi)
  test_gen(apt range)
  test_gen(apt bit_manip)

  # msperf_gen(apt block)
endif()
