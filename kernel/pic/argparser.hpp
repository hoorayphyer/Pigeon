#pragma once

#include <optional>

namespace pic {
struct CLIArgs {
  std::optional<std::string> journal_file{};
  std::optional<std::string> config_file{};
  std::optional<std::string> resume_dir{};
  bool is_dry_run = false;
};

auto parse_args(int argc, char** argv) {
  CLIArgs res;

  std::vector<std::string> args;
  {
    args.reserve(argc);
    for (int i = 0; i < argc; ++i) args.push_back({argv[i]});
  }

  for (int i = 1; i < args.size();) {
    if ("--dry-run" == args[i]) {
      res.is_dry_run = true;
      ++i;
    } else if ("--config" == args[i] or "-c" == args[i]) {
      res.config_file.emplace(args[i + 1]);
      i += 2;
    } else if ("--resume" == args[i]) {
      res.resume_dir.emplace(args[i + 1]);
      i += 2;
    } else if ("--journal" == args[i]) {
      res.journal_file.emplace(args[i + 1]);
      i += 2;
    }
  }
  return res;
}
}  // namespace pic
