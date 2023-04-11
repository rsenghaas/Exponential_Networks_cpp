#include <filesystem>
#include <iostream>
#include <string_view>
#include <vector>

#include "magic_numbers.h"

namespace fs = std::filesystem;

const std::string wsl_path = "/home/spamdoodler/SpamOS/Documents/Masterarbeit/Exponential_Networks_cpp";

const std::vector<std::string_view> directory_paths{
    "data/path_data", "data/map_data", "data/intersection_data", "graphics"};

auto create_directories() -> void {
  fs::path root_path = fs::path(getenv("HOME"));
  // std::string DataDirectory = fs::current_path();
  // root_path /= kDataDirectory;
  root_path = wsl_path;
  std::cout << fs::current_path() << std::endl;
  fs::current_path(root_path);  // (3)
  std::cout << "Current path is " << fs::current_path() << '\n';
  fs::remove_all("data");
  for (auto &p : directory_paths) {
    fs::create_directories(p);
  }
  // std::system("tree data");
}
