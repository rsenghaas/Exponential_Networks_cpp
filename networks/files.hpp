#include <filesystem>
#include <iostream>
#include <string_view>
#include <vector>

#include "magic_numbers.h"

namespace fs = std::filesystem;

const std::vector<std::string_view> directory_paths{
    "data/path_data", "data/map_data", "data/intersection_data", "graphics"};

auto create_directories() -> void {
  fs::path root_path = fs::path(getenv("HOME"));
  root_path /= kDataDirectory;
  fs::current_path(root_path);  // (3)
  std::cout << "Current path is " << fs::current_path() << '\n';
  for (auto &p : directory_paths) {
    fs::create_directories(p);
  }
  std::system("tree data");
}
