#ifndef _ARGPARSER_HPP_
#define _ARGPARSER_HPP_

namespace pic {
  struct CLIArgs {
    // std::string journal_file = "journal.txt";
    std::optional<std::string> parameters_file{};
    std::optional<std::string> checkpoint_dir {};
    bool is_screen = false;
  };

  auto parse_args( int argc, char** argv ) {
    CLIArgs res;

    std::vector<std::string> args;
    {
      args.reserve(argc);
      for ( int i = 0; i < argc; ++i )
        args.push_back({argv[i]});
    }

    for ( int i = 1; i < args.size(); ) {
      // if ( "--journal" == args[i] ) {
      //   res.journal_file = args[i + 1];
      //   i += 2;
      // } else
      if ( "--screen" == args[i] ) {
        res.is_screen = true;
        ++i;
        // ignore the rest of it
        break;
      }
      else if ( "--params" == args[i] ) {
        res.parameters_file.emplace( args[i + 1] );
        i += 2;
      }
      else if ( "--checkpoint" == args[i] ) {
        res.checkpoint_dir.emplace( args[i + 1] );
        i += 2;
      }
    }
    return res;
  }
}

#endif
