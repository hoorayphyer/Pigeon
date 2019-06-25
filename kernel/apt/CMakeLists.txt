if (${ENABLE_TESTS})
  test_gen(apt array) # DONE
  test_gen(apt pair) # DONE
  test_gen(apt "foreach") # DONE
  test_gen(apt block) # DONE
  test_gen(apt vec) # DONE
  # test_gen(apt virtual_vec) # TODO
  test_gen(apt numeric) # DONE
  test_gen(apt handle) # DONE
  test_gen(apt priority_queue) # DONE

  msperf_gen(apt block)
endif()
