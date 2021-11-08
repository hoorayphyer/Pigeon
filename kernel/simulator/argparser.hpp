#pragma once

#include <optional>
#include <string>
#include <vector>

namespace pic {
struct CLIArgs {
  std::optional<std::string> journal_file{};
  std::optional<std::string> config_file{};
  std::optional<std::string> resume_dir{};
  bool is_dry_run = false;
  std::vector<std::string> rest;
};

auto parse_args(int argc, char* argv[]) {
  CLIArgs res;

  // argc is 1 plus the number of arguments. But `argv[argc]` must be zero
  std::string val;
  for (int i = 0; i < argc;) {
    val = argv[i];
    if ("--dry-run" == val) {
      res.is_dry_run = true;
      ++i;
    } else if ("--config" == val or "-c" == val) {
      res.config_file.emplace(argv[i + 1]);
      i += 2;
    } else if ("--resume" == val) {
      res.resume_dir.emplace(argv[i + 1]);
      i += 2;
    } else if ("--journal" == val) {
      res.journal_file.emplace(argv[i + 1]);
      i += 2;
    } else {
      res.rest.push_back(std::move(val));
      ++i;
    }
  }

  return res;
}
}  // namespace pic
