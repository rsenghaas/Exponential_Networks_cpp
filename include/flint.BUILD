FLINT_HEADER_FILES = glob(
  ["**/*.h"] 
)

FLINT_SRC_FILES = glob(
  ["**/*.c"]
)

cc_library(
  name = "flint",
  hdrs = FLINT_HEADER_FILES,
  srcs = FLINT_SRC_FILES,
  includes = ["."],
  visibility = ["//visibility:public"],
)
