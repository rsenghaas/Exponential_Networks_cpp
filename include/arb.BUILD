ARB_HEADER_FILES = glob(
  ["**/*.h"] 
)

ARB_SRC_FILES = glob(
  ["**/*.c"]
)

cc_library(
  name = "arb",
  hdrs = ARB_HEADER_FILES,
  srcs = ARB_SRC_FILES,
  includes = ["."],
  visibility = ["//visibility:public"],
  deps = ["@flint_archive//:flint"],
)
